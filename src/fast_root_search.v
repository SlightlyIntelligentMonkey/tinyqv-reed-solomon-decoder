/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

//based off of https://doi.org/10.48550/arXiv.cs/0606035
module fast_root_search #(parameter MAX_ERRORS = 16)
                         (input clk, input rst, input [7:0] generator_polynomial, input [8*MAX_ERRORS-1:0] error_locator_flat,
                          input [7*8:0] reduction_matrix,
                          output wire done, output reg [8*MAX_ERRORS:0] roots);
    wire [7:0] error_locator [0:MAX_ERRORS-1];
    generate
        for (genvar j = 0; j < MAX_ERRORS-1; j = j + 1) begin : flatten_error_locator
            assign error_locator[j] = error_locator_flat[8*(j+1)-:8];
        end
    endgenerate
    
    localparam BITWIDTH = $clog2(MAX_ERRORS);
    localparam INNER_DEGREE = (MAX_ERRORS*2 - 4)%5 == 0 ? (MAX_ERRORS*2 - 4)/5 : (MAX_ERRORS*2 - 4)/5 + 1;

    reg [7:0] L [0:7];
    reg L_ready;
    reg [7:0] A [0:255];
    //reg [7:0] inner_degree;

    reg [7:0] i;
    reg [8:0] j;
    reg [7:0] alpha;
    reg setup_done;

    reg [BITWIDTH-1:0] root_counter;

    wire [7:0] gray_code_prev;
    reg [7:0] gray_code;
    wire [7:0] gray_code_bit;
    reg [2:0] delta;
    reg polynomial_started;
    reg [7:0] polynomial_result;

    //gray code is calculated as j ^ (j >> 1) then xor to find the bit that's different
    assign gray_code_prev = ((j[7:0]-1) ^ ((j[7:0]-1) >> 1));
    //assign gray_code_next = ((j[7:0]+1) ^ ((j[7:0]+1) >> 1));
    assign gray_code = (j[7:0] ^ (j[7:0] >> 1));
    assign gray_code_bit = gray_code ^ gray_code_prev;

    reg [7:0] mul_a [0:2];
    reg [7:0] mul_b [0:2];
    wire [7:0] mul_c [0:2];
    finite_field_multiplier_mastravito multiplier0(mul_a[0], mul_b[0], reduction_matrix, mul_c[0]);
    finite_field_multiplier_mastravito multiplier1(mul_a[1], mul_b[1], reduction_matrix, mul_c[1]);
    finite_field_multiplier_mastravito multiplier2(mul_a[2], mul_b[2], reduction_matrix, mul_c[2]);

    always @(posedge rst) begin
        i <= 0;
        j <= 0;
        root_counter <= 0;
        setup_done = 0;
        polynomial_result <= 0;
        polynomial_started <= 0;
        alpha <= generator_polynomial;
        //TODO: fix this later maybe
        //inner_degree = (MAX_ERRORS - 4)/5;
    end

    always @(posedge clk) begin
        //setup for fast evaluation search
        if(setup_done == 0 && j != 4) begin
            L[i[2:0]] <= L[i[2:0]] + error_locator[5*i + 2*j]*(alpha << j);
            j <= j + 1;
        end
        if(setup_done == 0 && j == 4 && i != INNER_DEGREE) begin
            mul_a[0] <= alpha;   mul_b[0] <= generator_polynomial;
            alpha <= mul_c[0];
            i <= i + 1;
            j <= 0;
        end
        if(setup_done == 0 && j == 4 && i == INNER_DEGREE) begin
            j <= 1;
            setup_done = 1;
        end
        if(setup_done == 1 && i == INNER_DEGREE && j != 256) begin

            integer k;
            for(k = 0; k < 8; k = k + 1) begin
                if (gray_code_bit[k] == 1)
                    delta = k[2:0];
            end

            A[j[7:0]] <= A[j[7:0]-1] ^ L[delta];

            //polynomial_result = error_locator[3]*gray_code*gray_code*gray_code;
            mul_a[0] <= gray_code;   mul_b[0] <= gray_code;
            mul_a[1] <= gray_code;   mul_b[1] <= error_locator[3];
            mul_a[2] <= mul_c[0];    mul_b[2] <= mul_c[1];
            polynomial_result <= mul_c[2];
            polynomial_started <= 1;

            j <= j + 1;
            i <= 0;
        end
        //evaluate inner polynomial
        if(setup_done == 1 && i != INNER_DEGREE) begin
            //temp = gray_code*gray_code;
            //temp = temp*temp;
            //temp = temp*gray_code;
            mul_a[0] <= gray_code;   mul_b[0] <= gray_code;
            mul_a[1] <= mul_c[0];    mul_b[1] <= mul_c[0];
            mul_a[2] <= mul_c[0];    mul_b[2] <= mul_c[1];

            polynomial_result <= polynomial_result ^ mul_c[2];
            i <= i + 1;
        end
        if(setup_done == 1 && polynomial_started == 1 && i == INNER_DEGREE) begin
            if (polynomial_result == 0) begin
                roots[root_counter] <= j[7:0];
                root_counter <= root_counter + 1;
            end
        end
    end
endmodule