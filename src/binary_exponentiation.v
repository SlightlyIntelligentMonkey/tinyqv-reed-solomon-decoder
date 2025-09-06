/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

module binary_exponentiation (input clk, input rst, input [7:0] element, input [7:0] power, input [7*8:0] reduction_matrix,
                              output data_ready, output reg [7:0] result);

    reg [7:0] exponent;
    reg [7:0] base;

    reg [7:0] mul_a;
    reg [7:0] mul_b;
    wire [7:0] mul_c;

    finite_field_multiplier_mastravito finite_field_multiplier(mul_a, mul_b, reduction_matrix, mul_c);

    assign data_ready = (exponent == 0) ? 1 : 0;
    always @(posedge rst) begin
        exponent = power;
        base <= element;
        result <= 1;
    end
    reg square_not_calc;
    always @(posedge clk) begin
        if (exponent != 0) begin
            if (square_not_calc == 0 && exponent[0] == 1) begin
                mul_a <= result;
                mul_b <= base;
                result <= mul_c;
                square_not_calc <= 1;
            end
            else begin
                mul_a <= base;  mul_b <= base;
                base <= mul_c;
                exponent = exponent >> 1;
            end
            if (square_not_calc == 1) begin
                square_not_calc <= 0;
            end
        end
    end
endmodule