
using BenchmarkExample

BenchmarkExample.PatchTest.𝐿 = 4.0

n = 16
# filename = "./msh/square_"*string(n)*".msh"
# BenchmarkExample.PatchTest.generateMsh(filename, transfinite=n+1)

filename = "./msh/square_irregular_"*string(n)*".msh"
BenchmarkExample.PatchTest.generateMsh(filename, lc = 0.25)