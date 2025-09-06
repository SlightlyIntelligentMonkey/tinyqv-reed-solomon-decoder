/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

//implementation of https://doi.org/10.1109/26.764911
module serial_berlekamp_massey #(parameter MAX_ERRORS = 16)
                                (input clk, input rst, input [$clog2(2*MAX_ERRORS)-1:0] code_length,
                                 input [8*(MAX_ERRORS*2)-1:0] syndromes_flat, input [7*8:0] reduction_matrix,
                                 output data_ready,
                                 output wire [8*MAX_ERRORS-1:0] out_error_locator,
                                 output wire [8*MAX_ERRORS-1:0] out_error_evaluator);

    wire [7:0] syndromes [0:(MAX_ERRORS*2)-1];
    generate
        for (genvar j = 0; j < (MAX_ERRORS*2)-1; j = j + 1) begin : unflatten_syndromes
            assign syndromes[j] = syndromes_flat[8*(j+1)-:8];
        end
    endgenerate

    localparam BITWIDTH = $clog2(2*MAX_ERRORS);
    reg [BITWIDTH-1:0] i;
    reg [BITWIDTH-1:0] j;
    reg error_locator_done;
    reg [7:0] auxillary_coeff;
    reg [7:0] auxillary_degree;

    reg [7:0] prev_discrepency;
    reg [7:0] discrepency;
    reg [7:0] delta;

    //1 is the current iteration 0 is the previous iteration
    reg [7:0] auxillary_polynomial [0:MAX_ERRORS-1];//[0:1];
    reg [7:0] error_locator [0:1][0:MAX_ERRORS-1];
    reg [7:0] error_evaluator [0:MAX_ERRORS-1];
    generate
        for (genvar j = 0; j < MAX_ERRORS-1; j = j + 1) begin : flatten_error_locator
            assign out_error_locator[8*(j+1)-:8] = error_locator[1][j];
        end
        for (genvar j = 0; j < MAX_ERRORS-1; j = j + 1) begin : flatten_error_evaluator
            assign out_error_evaluator[8*(j+1)-:8] = error_evaluator[j];
        end
    endgenerate

    reg [7:0] mul_a [0:2];
    reg [7:0] mul_b [0:2];
    wire [7:0] mul_c [0:2];
    finite_field_multiplier_mastravito multiplier0(mul_a[0], mul_b[0], reduction_matrix, mul_c[0]);
    finite_field_multiplier_mastravito multiplier1(mul_a[1], mul_b[1], reduction_matrix, mul_c[1]);
    finite_field_multiplier_mastravito multiplier2(mul_a[2], mul_b[2], reduction_matrix, mul_c[2]);

    integer k;
    always @(posedge rst) begin
        error_locator_done = 0;
        i <= 0;
        j <= 0;
        auxillary_degree <= 0;
        delta = 1;

        for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
            error_locator[0][k] = 1;
        end
        for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
            auxillary_polynomial[k] = 1;
        end

        discrepency <= syndromes[0];  //could be [1] too idk
    end

    always @(posedge clk) begin
        //calculating error locator polynomial
        if (j != code_length && error_locator_done == 0) begin
            //error_locator[1][j] = error_locator[j][0] * delta + discrepency * auxillary_polynomial[j - 1];
            //discrepency = discrepency + syndromes[i - j + 3]*error_locator[1][j - 1]
            mul_a[0] <= error_locator[0][j];     mul_b[0] <= delta;
            mul_a[1] <= discrepency;             mul_b[1] <= auxillary_polynomial[j - 1];
            mul_a[2] <= syndromes[i - j + 3];    mul_b[2] <= error_locator[1][j - 1];

            error_locator[1][j] <= mul_c[0] ^ mul_c[1];
            discrepency <= discrepency ^ mul_c[2];

            j <= j + 1;
        end
        //calculating error evaluator polynomial
        if (j != code_length && error_locator_done == 1) begin
            //error_evaluator[j] = error_evaluator[j - 1] + syndromes[i - j] * error_locator[1][j];
            mul_a[0] <= syndromes[i - j];    mul_b[0] <= error_locator[1][j];
            error_evaluator[j] <= error_evaluator[j - 1] ^ mul_c[0];

            j <= j + 1;
        end
        //control logic
        if (j == code_length && i != code_length) begin
            if (discrepency == 0 || auxillary_degree >= i + 1) begin
                auxillary_degree = auxillary_degree;
            end
            else begin
                auxillary_degree <= i + 1 - auxillary_degree;
                delta = discrepency;

                for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
                    auxillary_polynomial[k] = error_locator[0][k];
                end
            end
        
            if (error_locator_done == 0) begin
                //error_locator[0] = delta * error_locator[0][0];
                mul_a[0] = delta;   mul_b[0] = error_locator[0][0];
                error_locator[1][0] = mul_c[0];
                discrepency <= 0;
            end
            else begin
                //error_evaluator[0] = syndromes[i + 1] * error_locator[0][0];
                mul_a[1] <= syndromes[i + 1];   mul_b[0] <= error_locator[0][0];
                error_evaluator[0] <= mul_c[1];
            end

            //integer k;
            //for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
            //    error_locator[k][0] = error_locator[k][1];
            //end
            error_locator[0] = error_locator[1];//[0:MAX_ERRORS-1];
            //error_evaluator[0][0:MAX_ERRORS-1] = error_evaluator[1][0:MAX_ERRORS-1];

            j <= 1;
            i <= i + 1;
        end
        if (i == code_length && error_locator_done == 0) begin
            i <= 0;
            j <= 0;
            auxillary_degree <= 0;
            discrepency <= syndromes[0];  //could be [1] too idk
            for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
                auxillary_polynomial[k] = 1;
            end
            error_locator_done = 1;
        end
    end

    assign data_ready = (i == code_length && error_locator_done == 1) ? 1 : 0;

endmodule