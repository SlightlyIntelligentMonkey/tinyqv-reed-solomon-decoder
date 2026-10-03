/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */

module mastravito_reduction_matrix(input wire [8:0] irreducible_polynomial, output wire [7*8:0] reduction_matrix);
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
endmodule