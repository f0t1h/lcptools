#include "core.h"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

void log(const std::string &message) {
	std::cout << message << std::endl;
};

void test_core_constructors() {

	// test constructor with label, label_length, start, and end
	ublock *p1 = new ublock[1];
	p1[0] = 0b1111;
	struct core core1;
    init_core4(&core1, 4, p1, 2, 0, 10); // 1111 in binary

	assert(core1.bit_size == 4 && "Label length should be 4");
	assert(core1.bit_rep[0] == 0b1111 && "Label should be 0b1111");
	assert(core1.label == 2 && "Core label should be 1");
	assert(core1.start == 0 && "Start should be 0");
	assert(core1.end == 10 && "End should be 10");

    free_core(&core1);

	log("...  test_core_constructors passed!");
};

void test_core_compress() {

	ublock *p1 = new ublock[1];
	p1[0] = 0b101;
	ublock *p2 = new ublock[1];
	p2[0] = 0b111;
	struct core core1;
    init_core4(&core1, 3, p1, 10, 0, 3); // 101 in binary
	struct core core2;
    init_core4(&core2, 3, p2, 11, 0, 3); // 111 in binary
    
    core_compress(&core2, &core1);

	// expected result after compressing 101 and 111 is 10 (binary) => 2 in decimal
	assert(core1.bit_rep[0] == 0b10 && "Compressed core's label should be 0b10");
	assert(core1.bit_size == 2 && "Compressed core's label length should be 2");
	assert(core1.label == 10 && "Core's label should be 10");

    free_core(&core1);
    free_core(&core2);

	log("...  test_core_compress passed!");
};

void test_core_operator_overloads() {

	ublock *p1 = new ublock[1];
	p1[0] = 0b1010;
	ublock *p2 = new ublock[1];
	p2[0] = 0b1010;
	ublock *p3 = new ublock[1];
	p3[0] = 0b101;
	struct core core1;
    init_core4(&core1, 4, p1, 0, 0, 0); // 1010 in binary
	struct core core2;
    init_core4(&core2, 4, p2, 1, 1, 0); // 1010 in binary
	struct core core3;
    init_core4(&core3, 3, p3, 2, 2, 0); // 101 in binary

	assert(core_eq(&core1, &core2) && "core1 should be equal to core2");
	assert(core_neq(&core1, &core3) && "core1 should not be equal to core3");
	assert(core_lt(&core3, &core1) && "core3 should be less than core1");
	assert(core_gt(&core1, &core3) && "core1 should be greater than core3");
	assert(core_geq(&core1, &core2) && "core1 should be greater than or equal to core2");
	assert(core_leq(&core3, &core1) && "core3 should be less than or equal to core1");

    free_core(&core1);
    free_core(&core2);
    free_core(&core3);

	log("...  test_core_operator_overloads passed!");
};

int main() {
	log("Running test_core...");

	test_core_constructors();
	test_core_compress();
	test_core_operator_overloads();

	log("All tests in test_core completed successfully!");

	return 0;
};
