/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

module serial_syndrome_calculator #(parameter MAX_ERRORS = 16)
                            (input clk, input rst, input [7:0] generator_polynomial, input [8*256-1:0] coeff,
                             input [7*8:0] reduction_matrix,
                             output wire done, output wire [8*(MAX_ERRORS*2)-1:0] syndromes_flat);
    localparam COUNTERWIDTH = $clog2(2*MAX_ERRORS);
    reg evaluator_rst;
    reg [7:0] evaluator_x;
    wire evaluator_done;
    wire [7:0] evaluator_result;
    polynomial_evaluator evaluator(clk, evaluator_rst, evaluator_x, coeff, reduction_matrix, evaluator_done, evaluator_result);

    reg [7:0] mul_alpha;
    wire [7:0] mul_result;
    finite_field_multiplier_mastravito multiplier0(mul_alpha, generator_polynomial, reduction_matrix, mul_result);

    reg [7:0] syndromes [0:(MAX_ERRORS*2)-1];
    generate
        for (genvar j = 0; j < (MAX_ERRORS*2)-1; j = j + 1) begin : flatten_syndromes
            assign syndromes_flat[8*(j+1)-:8] = syndromes[j];
        end
    endgenerate

    reg [COUNTERWIDTH:0] counter;
    reg [7:0] alpha;
    assign done = (counter == MAX_ERRORS) ? 1 : 0;
    always @(posedge rst) begin
        counter <= 0;
        alpha <= generator_polynomial;
    end
    always @(posedge clk) begin
        evaluator_rst = 0;

        if (evaluator_done == 1 && counter != MAX_ERRORS) begin
            mul_alpha <= alpha;

            syndromes[counter[COUNTERWIDTH-1:0]] <= evaluator_result;
            evaluator_x <= mul_result;
            alpha <= mul_result;

            counter <= counter + 1;
        end
    end
    always @(negedge clk) begin
        if (evaluator_done == 1)
            evaluator_rst = 1;
    end
endmodule