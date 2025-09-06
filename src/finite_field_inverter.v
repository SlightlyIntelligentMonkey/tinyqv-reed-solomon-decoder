/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

//just use exponentiation
module finite_field_inverter (input clk, input rst, input [7:0] element, input [7*8:0] reduction_matrix,
                              output data_ready, output reg [7:0] inverse);

    binary_exponentiation power(clk, rst, element, 254, reduction_matrix, data_ready, inverse);
endmodule