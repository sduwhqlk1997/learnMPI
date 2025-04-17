#include "test.h"
#include <iostream>
test::test(/* args */)
{
    i=1;
}

test::~test()
{
}

void test::print_i(){
    std::cout<<i<<std::endl;
}