using Test, SumSpaces
using MathLink
@testset "mathematica_correction" begin
    uS = ([[0.]], [[0.]], [[0.]], [[0.]]); x = [0.0];
    (x1, x2, uS) = mathematica_corrections(1,1,1,x,uS,[1.0], 5, stabilise=true, xx1=[0.],xx2=[1.],maxrecursion=100)

end