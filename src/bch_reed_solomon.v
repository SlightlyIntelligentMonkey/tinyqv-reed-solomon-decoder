/*
 * Copyright (c) 2025 Dawson Hubbard
 * SPDX-License-Identifier: Apache-2.0
 */


//implementation of https://link.springer.com/chapter/10.1007/3-540-51083-4_67
//q is a reduction matrix
module finite_field_multiplier_mastravito (input [7:0] a, input [7:0] b, input [6:0] q [0:7], output wire [7:0] c);

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
                assign c_part[j][i+1] = d[i+8] & q[j][i];
            end
            assign c[j] = ^c_part[j];
        end
    endgenerate
endmodule

module finite_field_multiplier_serial ( input clk, input rst, input [7:0] a, input [7:0] b, input [8:0] q,
                                        output wire done, output reg [7:0] c);
    reg [8:0] bb;
    reg [4:0] counter;

    assign done = (counter == 8) ? 1 : 0;

    always @(posedge rst) begin
        c = 0;
        counter = 0;
        bb[7:0] = b;
    end
    always @(posedge clk) begin
        if (counter != 8) begin
            if (a[i] == 1)
                c = c ^ bb;
            bb = bb << 1;
            if (bb[8] & 1 == 1)
                bb = bb ^ q;
            counter = counter + 1;
        end
    end
endmodule

module binary_exponentiation (input clk, input rst, input [7:0] element, input [7:0] power, input [6:0] reduction_matrix [0:7],
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

//just use exponentiation
module finite_field_inverter (input clk, input rst, input [7:0] element, input [6:0] reduction_matrix [0:7],
                              output data_ready, output reg [7:0] inverse);

    binary_exponentiation power(clk, rst, element, 254, reduction_matrix, data_ready, inverse);
endmodule

//uses horners method
module polynomial_evaluator #(parameter DEGREE = 256)
                             (input clk, input rst, input [7:0] x, input [7:0] coeff [0:DEGREE-1], input [6:0] reduction_matrix [0:7],
                              output data_ready, output reg [7:0] result);
    localparam IWIDTH = $clog2(DEGREE);
    reg [IWIDTH:0] i;

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

module serial_syndrome_calculator #(parameter MAX_ERRORS = 16)
                            (input clk, input rst, input [7:0] generator_polynomial, input [7:0] coeff [0:255],
                             input [6:0] reduction_matrix [0:7],
                             output wire done, output reg [7:0] syndromes [0:(MAX_ERRORS*2)-1]);
    localparam COUNTERWIDTH = $clog2(2*MAX_ERRORS);
    reg evaluator_rst;
    reg [7:0] evaluator_x;
    wire evaluator_done;
    wire [7:0] evaluator_result;
    polynomial_evaluator evaluator(clk, evaluator_rst, evaluator_x, coeff, reduction_matrix, evaluator_done, evaluator_result);

    reg [7:0] mul_alpha;
    wire [7:0] mul_result;
    finite_field_multiplier_mastravito multiplier0(mul_alpha, generator_polynomial, reduction_matrix, mul_result);

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

//implementation of https://doi.org/10.1109/26.764911
module serial_berlekamp_massey #(parameter MAX_ERRORS = 16)
                                (input clk, input rst, input [$clog2(2*MAX_ERRORS)-1:0] code_length,
                                 input [7:0] syndrome [0:(MAX_ERRORS*2)-1], input [6:0] reduction_matrix [0:7],
                                 output data_ready,
                                 output wire [7:0] out_error_locator [0:MAX_ERRORS-1],
                                 output wire [7:0] out_error_evaluator [0:MAX_ERRORS-1]);
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

    assign out_error_locator = error_locator[1];
    assign out_error_evaluator = error_evaluator;

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

        discrepency <= syndrome[0];  //could be [1] too idk
    end

    always @(posedge clk) begin
        //calculating error locator polynomial
        if (j != code_length && error_locator_done == 0) begin
            //error_locator[1][j] = error_locator[j][0] * delta + discrepency * auxillary_polynomial[j - 1];
            //discrepency = discrepency + syndrome[i - j + 3]*error_locator[1][j - 1]
            mul_a[0] <= error_locator[0][j];     mul_b[0] <= delta;
            mul_a[1] <= discrepency;             mul_b[1] <= auxillary_polynomial[j - 1];
            mul_a[2] <= syndrome[i - j + 3];     mul_b[2] <= error_locator[1][j - 1];

            error_locator[1][j] <= mul_c[0] ^ mul_c[1];
            discrepency <= discrepency ^ mul_c[2];

            j <= j + 1;
        end
        //calculating error evaluator polynomial
        if (j != code_length && error_locator_done == 1) begin
            //error_evaluator[j] = error_evaluator[j - 1] + syndrome[i - j] * error_locator[1][j];
            mul_a[0] <= syndrome[i - j];     mul_b[0] <= error_locator[1][j];
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
                //error_evaluator[0] = syndrome[i + 1] * error_locator[0][0];
                mul_a[1] <= syndrome[i + 1];   mul_b[0] <= error_locator[0][0];
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
            discrepency <= syndrome[0];  //could be [1] too idk
            for (k = 0; k < MAX_ERRORS-1; k = k + 1) begin
                auxillary_polynomial[k] = 1;
            end
            error_locator_done = 1;
        end
    end

    assign data_ready = (i == code_length && error_locator_done == 1) ? 1 : 0;

endmodule

//based off of https://doi.org/10.48550/arXiv.cs/0606035
module fast_root_search #(parameter MAX_ERRORS = 16)
                         (input clk, input rst, input [7:0] generator_polynomial, input [7:0] error_locator [0:MAX_ERRORS-1],
                          input [6:0] reduction_matrix [0:7],
                          output wire done, output reg [7:0] roots [0:MAX_ERRORS-1]);
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

module forney_algorithm #(parameter MAX_ERRORS = 16)
                         (input clk, input rst, input [7:0] first_root, input [7:0] roots[0:MAX_ERRORS-1],
                          input [7:0] error_locator [0:MAX_ERRORS-1], input [7:0] error_evaluator [0:MAX_ERRORS-1],
                          input [6:0] reduction_matrix [0:7],
                          output done, input [7:0] message_data [0:255], output reg [7:0] corrected_data [0:255]);
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
        evaluator_evaluator(clk, evaluator_rst, evaluator_x, error_evaluator, reduction_matrix, evaluator_done, evaluator_result);

    reg [7:0] error_locator_derivative [0:MAX_ERRORS-1];
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

        for (k = 0; k < 256; k = k + 1) begin
            corrected_data[k] = message_data[k];
        end
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
                corrected_data[error_position[2]] <= corrected_data[error_position[2]] ^ mul_c;

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