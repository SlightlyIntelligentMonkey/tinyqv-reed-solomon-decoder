/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */


//implementation of https://link.springer.com/chapter/10.1007/3-540-51083-4_67
//q is a reduction matrix
module finite_field_multiplier_mastravito (input [7:0] a, input [7:0] b, input [7*8:0] q, output wire [7:0] c);

    wire [7:0] d_part [0:14];
    wire [14:0] d;
    wire [7:0] c_part [0:7];

    //multiplication
    generate
        for (genvar i = 0; i < 15; i = i + 1) begin : multiplication
            for (genvar j = 0; j < 8; j = j + 1) begin
                if (i - j > 0 && i - j < 8)
                    assign d_part[i][j] = a[i - j] & b[j];
                else
                    assign d_part[i][j] = 0;
            end
            assign d[i] = ^d_part[i];
        end

        //reduction
        for (genvar j = 0; j < 8; j = j + 1) begin : reduction_init
            assign c_part[j][0] = d[j];
        end
        for (genvar j = 0; j < 8; j = j + 1) begin : reduction
            for (genvar i = 0; i < 7; i = i + 1) begin
                //assign c_part[j][i+1] = d[i+8] & q[j][i];
                assign c_part[j][i+1] = d[i+8] & q[7*j+i];
            end
            assign c[j] = ^c_part[j];
        end
    endgenerate
endmodule

//serial alternative
module finite_field_multiplier_serial ( input clk, input rst, input [7:0] a, input [7:0] b, input [8:0] q,
                                        output wire done, output reg [7:0] c);
    reg [8:0] bb;
    reg [4:0] i;

    assign done = (i == 8) ? 1 : 0;

    always @(posedge rst) begin
        c = 0;
        i = 0;
        bb[7:0] = b;
    end
    always @(posedge clk) begin
        if (i != 8) begin
            if (a[i] == 1)
                c = c ^ bb;
            bb = bb << 1;
            if (bb[8] & 1 == 1)
                bb = bb ^ q;
            i = i + 1;
        end
    end
endmodule