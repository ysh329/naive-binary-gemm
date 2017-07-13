#!/bin/python

def popcnt(n1, n2):
    int2bin_list = lambda num: map(int, bin(num).replace("0b", ""))

    n1_list = int2bin_list(n1)
    n2_list = int2bin_list(n2)

    if len(n1_list) < len(n2_list):
        zero_list = [0] * (len(n2_list)-len(n1_list))
        zero_list.extend(n1_list)
        n1_list = zero_list
    elif len(n1_list) > len(n2_list):
        zeros_list = [0] * (len(n1_list)-len(n2_list)).extend(n2_list)
        zero_list.extend(n2_list)
        n2_list = zero_list
    #########
    # Check #
    #########
    print n1_list
    print n2_list

    xor_res_list = map(lambda nn1, nn2: nn1^nn2, n1_list, n2_list)

    return xor_res_list.count(1)

if __name__ == "__main__":
    print popcnt(1,1)
    print popcnt(1,2)
