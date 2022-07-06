using Test, SumSpacesMathLink

@testset "mathematica_correction" begin
    
    λ=1; μ=1; η=1

    if isdir("uS-lmbda-$λ-mu-$μ-eta-$η")
        rm("uS-lmbda-$λ-mu-$μ-eta-$η", recursive=true)
    end

    uS = ([[0.]], [[0.]], [[0.]], [[0.]]); x = [0.0];
    (x1, x2, uS) = mathematica_corrections(λ, μ, η, x, uS, [1.0], 5, stabilise=true, xx1=[0.], xx2=[1.], maxrecursion=100)

    @test x1 == [0.0]
    @test x2 == [0.0, 1.0]
    @test uS == ([[0.6470543404462853]], [[-0.22602990220463598]], [[0.0, 0.07303314441297332]], [[0.0, 0.07303129615272594]])

    rm("uS-lmbda-$λ-mu-$μ-eta-$η", recursive=true)

    uS = ([[0.]], [[0.]], [[0.]], [[0.]]); x = [0.0];
    (x1, x2, uS) = mathematica_corrections(λ, μ, η, x, uS, [1.0], 5, stabilise=false, xx1=[0.], xx2=[1.], maxrecursion=100)
    @test x1 == [0.0]
    @test x2 == [0.0, 1.0]
    @test uS == ([[0.6470543404462853]], [[-0.22602990220463598]], [[0.0, 0.7158702689927211]], [[0.0, 0.06270043552531117]])

    rm("uS-lmbda-$λ-mu-$μ-eta-$η", recursive=true)

    uS = ([[0.1]], [[0.2]], [[0.3]], [[0.4]]); x = [1.0];
    (x1, x2, uS1) = mathematica_corrections(λ, μ, η, x, uS,[1.0], 5, stabilise=true, xx1=[0.], xx2=[1.], maxrecursion=100)
    @test x1 == [0.0, 1.0]
    @test x2 == [1.0]
    @test uS1 == ([[0.6470543404462853, 0.1]], [[-0.22602990220463598, 0.2]], [[0.07303314441297332]], [[0.07303129615272594]])

    
    (x3, x4, uS2) = mathematica_corrections(λ, μ, η, x, uS ,[1.0], 5, stabilise=true, xx1=[0.], xx2=[1.], maxrecursion=100)
    @test x3 == x1
    @test x4 == x2
    @test uS2 == uS1

    rm("uS-lmbda-$λ-mu-$μ-eta-$η", recursive=true)

    uS = ([0.1.*ones(5)], [0.2.*ones(5)], [0.3.*ones(5)], [0.4.*ones(5)]); x = 1.0:1.0:5;
    (x1, x2, uS1) = mathematica_corrections(λ, μ, η, x, uS, [1.0], 5, stabilise=true, xx1=[0.], xx2=[1.], maxrecursion=100)
    @test x1 == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    @test x2 == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test uS1 == ([[0.6470543404462853, 0.1, 0.1, 0.1, 0.1, 0.1]], [[-0.22602990220463598, 0.2, 0.2, 0.2, 0.2, 0.2]], [[0.07303314441297332, 0.3, 0.3, 0.3, 0.3]], [[0.07303129615272594, 0.4, 0.4, 0.4, 0.4]])
    
    # Check the saved-load trick works
    (x1, x2, uS2) = mathematica_corrections(λ, μ, η, x, uS, [1.0], 5, stabilise=true, xx1=[0.], xx2=[1.], maxrecursion=100)
    @test x1 == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    @test x2 == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test uS2 == uS1

    rm("uS-lmbda-$λ-mu-$μ-eta-$η", recursive=true)


end