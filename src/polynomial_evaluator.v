/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

//uses horners method
module polynomial_evaluator #(parameter DEGREE = 256)
                             (input clk, input rst, input [7:0] x, input [8*DEGREE-1:0] coeff_flat, input [7*8:0] reduction_matrix,
                              output data_ready, output reg [7:0] result);
    localparam IWIDTH = $clog2(DEGREE);
    reg [IWIDTH:0] i;

    wire [7:0] coeff [0:DEGREE-1];
    generate
        for (genvar j = 0; j < DEGREE-1; j = j + 1) begin : unflatten_coeff
            assign coeff[j] = coeff_flat[8*(j+1)-:8];
        end
    endgenerate

    assign data_ready = (i == ~0) ? 1 : 0;

    always @(posedge rst) begin
        i <= DEGREE - 1;
        result <= 0;
    end

    reg [7:0] mul_result;
    finite_field_multiplier_mastravito multiplier(x, result, reduction_matrix, mul_result);
    always @(posedge clk) begin
        if (i != ~0) begin
            result <= mul_result ^ coeff[i[IWIDTH-1:0]];
            //result = result + result*x + coeff[i]; 
            i <= i - 1;
        end
    end
endmodule