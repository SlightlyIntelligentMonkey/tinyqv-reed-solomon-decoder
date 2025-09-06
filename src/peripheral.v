/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

`default_nettype none

module tqvp_reed_solomon_decoder (
    input         clk,          // Clock - the TinyQV project clock is normally set to 64MHz.
    input         rst_n,        // Reset_n - low to reset.

    input  [7:0]  ui_in,        // The input PMOD, always available.  Note that ui_in[7] is normally used for UART RX.
                                // The inputs are synchronized to the clock, note this will introduce 2 cycles of delay on the inputs.

    output [7:0]  uo_out,       // The output PMOD.  Each wire is only connected if this peripheral is selected.
                                // Note that uo_out[0] is normally used for UART TX.

    input [5:0]   address,      // Address within this peripheral's address space
    input [31:0]  data_in,      // Data in to the peripheral, bottom 8, 16 or all 32 bits are valid on write.

    // Data read && write requests from the TinyQV core.
    input [1:0]   data_write_n, // 11 = no write, 00 = 8-bits, 01 = 16-bits, 10 = 32-bits
    input [1:0]   data_read_n,  // 11 = no read,  00 = 8-bits, 01 = 16-bits, 10 = 32-bits
    
    output [31:0] data_out,     // Data out from the peripheral, bottom 8, 16 or all 32 bits are valid on read when data_ready is high.
    output        data_ready,

    output        user_interrupt  // Dedicated interrupt request for this peripheral
);
    reg [7:0] block_length;
    reg [7:0] message_length;
    reg [7:0] code_length;

    reg [7:0] generator_polynomial;
    reg [8:0] irreducible_polynomial;
    reg [7:0] first_root;

    //calculate reduction matrix
    /* verilator lint_off UNOPTFLAT */
    reg [7*8:0] reduction_matrix;
    generate
        genvar j;
        genvar i;
        for (j = 0; j < 8; j = j + 1)
            assign reduction_matrix[7*j] = irreducible_polynomial[j];

        for (j = 0; j < 8; j = j + 1) begin
            for (i = 1; i < 7; i = i + 1) begin
                if (j - 1 >= 0)
                    assign reduction_matrix[7*j + i] = reduction_matrix[7*(j - 1) + i - 1] ^ reduction_matrix[7*7 + i - 1];
                else
                    assign reduction_matrix[7*j + i] = reduction_matrix[7*7 + i - 1];
            end
        end        
    endgenerate
    /* verilator lint_on UNOPTFLAT */


    reg [7:0] accum;

    //32 bit words addressed with address range 0 - 63
    reg [8*256-1:0] message_data [0:2];
    wire [8*256-1:0] decoded_data;
    wire [8*256-1:0] error_polynomial;

    //decode reed solomon code
    localparam MAX_ERRORS = 16;
    reg [8*(MAX_ERRORS*2)-1:0] syndromes;

    wire syndrome_done;
    wire berlekamp_massey_done;
    wire root_search_done;
    wire forney_algorithm_done;

    wire syndrome_rst;
    wire berlekamp_massey_rst;
    wire root_search_rst;
    wire forney_algorithm_rst;

    wire [8*(MAX_ERRORS*2)-1:0] calculated_syndromes;
    serial_syndrome_calculator #(MAX_ERRORS)
        syndrome_calculator(clk, syndrome_rst, generator_polynomial, message_data[0],reduction_matrix,
        syndrome_done, calculated_syndromes);

    wire [7:0] berlekamp_massey_code_length;
    reg [8*MAX_ERRORS-1:0] error_locator [0:2];
    reg [8*MAX_ERRORS-1:0] error_evaluator [0:2];
    serial_berlekamp_massey #(MAX_ERRORS)
        berlekamp_massey(clk, berlekamp_massey_rst, berlekamp_massey_code_length, syndromes, reduction_matrix,
                         berlekamp_massey_done, error_locator[0], error_evaluator[0]);
    
    wire [8*MAX_ERRORS-1:0] root_search_roots;
    fast_root_search #(MAX_ERRORS)
        root_search(clk, root_search_rst, generator_polynomial, error_locator[1], reduction_matrix,
                    root_search_done, root_search_roots);

    forney_algorithm #(MAX_ERRORS)
        forney(clk, forney_algorithm_rst, first_root, root_search_roots, error_locator[2], error_evaluator[2], reduction_matrix,
               forney_algorithm_done, message_data[2], decoded_data);

    assign data_ready = (data_read_n == 'b10) ? 1 : 0;
    always @(posedge clk) begin

        if (data_write_n == 'b10 && ui_in[0] == 0) begin
            message_data[0][8*((address << 2) + 1)-:8] = data_in[7:0];
            message_data[0][8*((address << 2) + 2)-:8] = data_in[15:8];
            message_data[0][8*((address << 2) + 3)-:8] = data_in[23:16];
            message_data[0][8*((address << 2) + 4)-:8] = data_in[31:24];
        end
        if (data_write_n != 'b11 && address == 0 && ui_in[0] == 1) begin
            block_length = data_in[7:0];
            code_length = block_length - message_length;
            $display("Wrote to block_length: %d", block_length);
        end
        if (data_write_n != 'b11 && address == 1 && ui_in[0] == 1) begin
            message_length = data_in[7:0];
            code_length = block_length - message_length;
            $display("Wrote to message_length: %d", message_length);
        end
        if (data_write_n != 'b11 && address == 2 && ui_in[0] == 1) begin
            generator_polynomial = data_in[7:0];
            $display("Wrote to generator_polynomial: %d", generator_polynomial);
        end
        if (data_write_n == 'b01 && address == 3 && ui_in[0] == 1) begin
            irreducible_polynomial = data_in[8:0];
            $display("Wrote to irreducible_polynomial: %d", irreducible_polynomial);
        end
        if (data_write_n != 'b11 && address == 4 && ui_in[0] == 1) begin
            first_root = data_in[7:0];
            $display("Wrote to first_root: %d", first_root);
        end

        //if (data_read_n == 'b11)
        //    data_ready = 0;
        if(data_read_n == 'b10) begin
            data_out[7:0]   = decoded_data[8*((address << 2) + 1)-:8];
            data_out[15:8]  = decoded_data[8*((address << 2) + 2)-:8];
            data_out[23:16] = decoded_data[8*((address << 2) + 3)-:8];
            data_out[31:24] = decoded_data[8*((address << 2) + 4)-:8];
        //    data_ready = 1;
        end

        //syndrome_rst <= 0;
        //berlekamp_massey_rst <= 0;
        //root_search_rst <= 0;
        //forney_algorithm_rst <= 0;

        //go to next pipeline stage
        if (ui_in[1] == 1 && syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) begin
            message_data[2] <= message_data[1];
            message_data[1] <= message_data[0];
            syndromes <= calculated_syndromes;

            error_locator[2] <= error_locator[1];
            error_locator[1] <= error_locator[0];
            error_evaluator[2] <= error_evaluator[1];
            error_evaluator[1] <= error_evaluator[0];

            //syndrome_rst <= 1;
            //berlekamp_massey_rst <= 1;
            //root_search_rst <= 1;
            //forney_algorithm_rst <= 1;
        end
    end

    always @(negedge rst_n) begin
        $display("RESET!!!!");
        block_length = 255;
        message_length = 223;
        code_length = 32;

        first_root = 41;
        generator_polynomial = 79;
        irreducible_polynomial = 'b111110101;
    end

    assign syndrome_rst = (ui_in[1] == 1 && syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) ||rst_n == 0 ? 1 : 0;
    assign berlekamp_massey_rst = (ui_in[1] == 1 && syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) || rst_n == 0 ? 1 : 0;
    assign root_search_rst = (ui_in[1] == 1 && syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) || rst_n == 0 ? 1 : 0;
    assign forney_algorithm_rst = (ui_in[1] == 1 && syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) || rst_n == 0 ? 1 : 0;
    assign uo_out[0] = syndrome_done;
    assign uo_out[1] = berlekamp_massey_done;
    assign uo_out[2] = root_search_done;
    assign uo_out[3] = forney_algorithm_done;
    //assign user_interrupt = (syndrome_done == 1 && berlekamp_massey_done == 1 && root_search_done == 1 && forney_algorithm_done == 1) ? 1 : 0;
endmodule
