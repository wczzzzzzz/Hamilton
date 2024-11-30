
using BenchmarkExample.PatchTest

PatchTest.𝐿 = 4.0

n = 10
filename = "./msh/square_"*string(n)*".msh"

PatchTest.generateMsh(filename, transfinite=n+1)