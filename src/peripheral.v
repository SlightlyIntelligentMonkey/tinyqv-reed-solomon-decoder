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

    // Data read and write requests from the TinyQV core.
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
    reg [6:0] reduction_matrix [7:0];
    reg [7:0] first_root;

    reg [7:0] accum;

    //32 bit words addressed with address range 0 - 63
    reg [7:0] message_data [0:255];
    reg [7:0] decoded_data [0:255];

    //decode reed solomon code
    localparam MAX_CODE_LENGTH = 32;
    reg [7:0] syndromes [0:MAX_CODE_LENGTH];
    reg [7:0] located_data [0:255];
    reg [7:0] corrected_data [0:255];

    wire syndrome_done;
    wire berlekamp_massey_done;
    wire root_search_done;
    wire forney_algorithm_done;

    wire syndrome_rst;
    wire berlekamp_massey_rst;
    wire root_search_rst;
    wire forney_algorithm_rst;

    wire [7:0] calculated_syndromes [0:MAX_CODE_LENGTH]
    serial_syndrome_calculator #(MAX_CODE_LENGTH = MAX_CODE_LENGTH)
        syndrome_calculator(clk, syndrome_rst, generator_polynomial, message_data, reduction_matrix,
        syndrome_done, calculated_syndromes);

    wire [7:0] berlekamp_massey_code_length;
    wire [7:0] berlekamp_massey_error_locator [0:MAX_CODE_LENGTH/2];
    wire [7:0] berlekamp_massey_error_evaluator [0:MAX_CODE_LENGTH/2];
    serial_berlekamp_massey #(MAX_CODE_LENGTH = MAX_CODE_LENGTH)
        berlekamp_massey(clk, berlekamp_massey_rst, berlekamp_massey_code_length, syndromes, reduction_matrix,
                         berlekamp_massey_done, berlekamp_massey_error_locator, berlekamp_massey_error_evaluator);
    
    wire [7:0] root_search_error_locator [0:MAX_CODE_LENGTH/2];
    wire [7:0] root_search_roots [0:MAX_CODE_LENGTH];
    fast_root_search #(MAX_CODE_LENGTH = MAX_CODE_LENGTH)
        root_search(clk, root_search_rst, generator_polynomial, root_search_error_locator, reduction_matrix,
                    root_search_done, root_search_roots);

    wire [7:0] forney_error_locator_derivative [0:MAX_CODE_LENGTH];
    wire [7:0] forney_error_evaluator [0:MAX_CODE_LENGTH];
    forney_algorithm #(MAX_CODE_LENGTH = MAX_CODE_LENGTH)
        forney(clk, forney_algorithm_rst, first_root, forney_error_locator_derivative, forney_error_evaluator, reduction_matrix,
               forney_algorithm_done, corrected_data);

    always @(posedge clk) begin

        if (data_write_n == b10 and ui_in[0] == 0) begin
            message_data[(address << 2) + 0] = data_in[0:7];
            message_data[(address << 2) + 1] = data_in[7:15];
            message_data[(address << 2) + 2] = data_in[15:23];
            message_data[(address << 2) + 3] = data_in[23:31];
        end
        if (data_write_n != b11 and address == 0 and ui_in[0] == 1) begin
            block_length = data_in[0:7];
            code_length = block_length - message_length;
        end
        if (data_write_n != b11 and address == 1 and ui_in[0] == 1) begin
            message_length = data_in[0:7];
            code_length = block_length - message_length;
        end
        if (data_write_n != b11 and address == 2 and ui_in[0] == 1) begin
            generator_polynomial = data_in[0:7];
        end
        if (data_write_n == b01 and address == 3 and ui_in[0] == 1) begin
            irreducible_polynomial = data_in[0:8];
            //calculate reduction matrix
        end
        if (data_write_n != b11 and address == 4 and ui_in[0] == 1) begin
            first_root = data_in[0:7];
        end

        if (data_read_n == b11)
            data_ready = 0;
        if(data_read_n == b10) begin
            data_out[0:7]   = decoded_data[(address << 2) + 0];
            data_out[7:15]  = decoded_data[(address << 2) + 1];
            data_out[15:23] = decoded_data[(address << 2) + 2];
            data_out[23:31] = decoded_data[(address << 2) + 3];
            data_ready = 1;
        end

        syndrome_rst <= 0;
        berlekamp_massey_rst <= 0;
        root_search_rst <= 0;
        forney_algorithm_rst <= 0;

        //go to next pipeline stage
        if (syndrome_done == 1 and berlekamp_massey_done == 1 and root_search_done == 1 and forney_algorithm_done == 1) begin
            decoded_data <= corrected_data;
            corrected_data <= located_data;
            located_data <= message_data;
            syndromes <= calculated_syndromes;

            syndrome_rst <= 1;
            berlekamp_massey_rst <= 1;
            root_search_rst <= 1;
            forney_algorithm_rst <= 1;
        end
    end
endmodule
