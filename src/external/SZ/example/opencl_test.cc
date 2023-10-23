#include <sz.h>
#include <iostream>

int main(int argc, char *argv[])
{
	struct sz_opencl_state* state = nullptr;

	int error = sz_opencl_init(&state);
	if(error) {
		const char* msg = sz_opencl_error_msg(state);
		int code = sz_opencl_error_code(state);
		std::cerr << "ERROR: " << msg << " " << code << std::endl;
		goto done;
	}

	error |= sz_opencl_check(state);
	if(error) {
		const char* msg = sz_opencl_error_msg(state);
		int code = sz_opencl_error_code(state);
		std::cerr << "ERROR: " << msg << " " << code << std::endl;
	}

done:
	error |= sz_opencl_release(&state);

	std::cout << ((error)?"failed":"passed") << std::endl;
	return error;
}
