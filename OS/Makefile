obj = VRR_code_generator.o # vrr_code_generator.o

vrrcg: $(obj)
	g++ -std=c++17 -o vrrcg $(obj)

VRR_code_generator.o : VRR_code_generator.cxx
	g++ -std=c++17 -c VRR_code_generator.cxx
# vrr_code_generator.o: vrr_code_generator.cxx 
.PHONY : clean
clean :
	rm vrrcg $(obj)

