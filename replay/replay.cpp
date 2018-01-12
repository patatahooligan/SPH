// replay.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include "loadstate.h"

static LoadState *load_state_pointer;

void idle_func() {
	assert(load_state_pointer);

}

int main(int argc, char **argv) {
	if (argc != 2)
		throw std::runtime_error("Wrong number of arguments");
	LoadState load_state(argv[1]);
	load_state_pointer = &load_state;
    return 0;
}