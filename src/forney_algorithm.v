/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

module forney_algorithm #(parameter MAX_ERRORS = 16)
                         (input clk, input rst, input [7:0] first_root, input [8*MAX_ERRORS-1:0] roots_flat,
                          input [8*MAX_ERRORS-1:0] error_locator_flat, input [8*MAX_ERRORS-1:0] error_evaluator_flat,
                          input [7*8:0] reduction_matrix,
                          output done, input [8*256-1:0] message_data_flat, output reg [8*256-1:0] corrected_data_flat);
    wire [7:0] roots [0:MAX_ERRORS-1];
    generate
        for (genvar j = 0; j < MAX_ERRORS-1; j = j + 1) begin : unflatten_roots_flat
            assign roots[j] = roots_flat[8*(j+1)-:8];
        end
    endgenerate

    localparam BITWIDTH = $clog2(MAX_ERRORS);

    reg ffi_rst [0:1];
    wire ffi_done [0:1];
    reg [7:0] ffi_input [0:1];
    wire [7:0] ffi_result [0:1];

    finite_field_inverter ffi0(clk, ffi_rst[0], ffi_input[0], reduction_matrix, ffi_done[0], ffi_result[0]);
    finite_field_inverter ffi1(clk, ffi_rst[1], ffi_input[1], reduction_matrix, ffi_done[1], ffi_result[1]);

    reg power_rst;
    wire power_done;
    reg [7:0] power_input;
    reg [7:0] power_exponent;
    wire [7:0] power_result;
    binary_exponentiation power(clk, power_rst, power_input, power_exponent, reduction_matrix, power_done, power_result);

    reg [7:0] mul_a;
    reg [7:0] mul_b;
    wire [7:0] mul_c;
    finite_field_multiplier_mastravito multiplier0(mul_a, mul_b, reduction_matrix, mul_c);

    reg evaluator_rst;
    reg [7:0] evaluator_x;
    //wire [7:0] evaluator_coeff;
    wire evaluator_done;
    wire [7:0] evaluator_result;
    polynomial_evaluator #(MAX_ERRORS)
        evaluator_evaluator(clk, evaluator_rst, evaluator_x, error_evaluator_flat, reduction_matrix, evaluator_done, evaluator_result);

    //reg [7:0] error_locator_derivative [0:MAX_ERRORS-1];
    reg [8*MAX_ERRORS-1:0] error_locator_derivative;
    reg locator_rst;
    reg [7:0] locator_x;
    //wire [7:0] locator_coeff;
    wire locator_done;
    wire [7:0] locator_result;
    polynomial_evaluator #(MAX_ERRORS)
        locator_evaluator(clk, locator_rst, locator_x, error_locator_derivative, reduction_matrix, locator_done, locator_result);
    
    reg [7:0] error_position [0:2];
    reg [7:0] to_the_power [0:2];
    reg [7:0] locator_result_prev;
    //reg [7:0] ffi_locator_result [0:1];

    reg derivative_done;
    reg prev_j_lsb;
    reg [BITWIDTH:0] j;
    integer k;
    always @(posedge rst) begin
        j <= 0;
        power_exponent = 1 - first_root;
        derivative_done = 0;

        corrected_data_flat <= message_data_flat;
    end

    always @(posedge clk) begin
        //reset resets
        ffi_rst[0] <= 0;
        ffi_rst[1] <= 0;
        power_rst <= 0;
        locator_rst <= 0;
        evaluator_rst <= 0;

        prev_j_lsb = j[0];

        //calculate derivative of error locator
        if (derivative_done == 0 && j != MAX_ERRORS) begin
            error_locator_derivative[j[BITWIDTH-1:0]] <= 0;
            j <= j + 1;
        end
        if (derivative_done == 0 && j == MAX_ERRORS) begin
            derivative_done = 1;
            j <= 0;
        end
        if (derivative_done == 1 && j == 0) begin
            ffi_input[0] <= roots[0];
            power_input <= roots[0];

            ffi_rst[0] <= 1;
            power_rst <= 1;

            j <= j + 1;
        end
        //correct message data
        if (derivative_done == 1 && j != 0 && j != MAX_ERRORS) begin
            if (prev_j_lsb != j[0]) begin
                ffi_input[0] <= roots[j[BITWIDTH-1:0]];
                power_input <= roots[j[BITWIDTH-1:0]];
                ffi_rst[0] <= 1;
                power_rst <= 1;

                evaluator_x <= error_position[0];
                locator_x <= error_position[0];
                evaluator_rst <= 1;
                locator_rst <= 1;

                ffi_input[1] <= locator_result_prev;
                ffi_rst[1] <= 1;
            end

            if(ffi_done[0] == 1 && ffi_done[1] == 1 && evaluator_done == 1 && locator_done == 1 && power_done == 1) begin
                //correction = power_result * evaluator_result * ffi_result[1]
                mul_a <= evaluator_result;  mul_b <= ffi_result[1];
                mul_a <= mul_c;             mul_b <= to_the_power[2];
                corrected_data_flat[8*(error_position[2]+1)-:8] <= corrected_data_flat[8*(error_position[2]+1)-:8] ^ mul_c;

                error_position[2] <= error_position[1];
                error_position[1] <= error_position[0];
                error_position[0] <= ffi_result[0];

                to_the_power[2] <= to_the_power[1];
                to_the_power[1] <= to_the_power[0];
                to_the_power[0] <= power_result;

                locator_result_prev <= locator_result;
                j <= j + 1;
            end
        end
    end
    assign done = (derivative_done == 1 && j == MAX_ERRORS) ? 1 : 0;
endmodule