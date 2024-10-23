# Compile & run:
compile: g++ -o some_name cpp_files.cpp 

run: ./some_name  path_to_your_graph_file.txt

# To extend the stack:
ulimit -s unlimited 
