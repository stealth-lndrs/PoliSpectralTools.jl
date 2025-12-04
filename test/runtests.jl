using Test
using LinearAlgebra
using SparseArrays
using SpectralTools

@testset "SpectralTools" begin
    @testset "Chebyshev nodes" begin
        x = cheb_lobatto_nodes(10)
        @test length(x) == 11
        @test isapprox(x[1], 1.0; atol=eps())
        @test isapprox(x[end], -1.0; atol=eps())
        @test all(diff(x) .< 0)
    end

    @testset "Differentiation matrices" begin
        D, D2 = cheb_D_matrices(12)
        x = cheb_lobatto_nodes(12)
        f = x .^ 3
        df = 3 .* x .^ 2
        d2f = 6 .* x
        @test maximum(abs.(D * f - df)) < 1e-10
        @test maximum(abs.(D2 * f - d2f)) < 1e-10
    end

    @testset "Laplacian operator" begin
        _, D2x = cheb_D_matrices(5)
        _, D2y = cheb_D_matrices(7)
        A = laplacian_plus_I_operator(D2x, D2y)
        nx = size(D2x, 1)
        ny = size(D2y, 1)
        @test size(A, 1) == nx * ny
        @test size(A, 2) == nx * ny
        @test issparse(A)
        ix = 2
        iy = 3
        idx = ix + (iy - 1) * nx
        @test isapprox(A[idx, idx], D2x[ix, ix] + D2y[iy, iy] + 1.0; atol=1e-10)
    end

    @testset "Dirichlet enforcement" begin
        nx = 5
        ny = 4
        _, D2x = cheb_D_matrices(nx - 1)
        _, D2y = cheb_D_matrices(ny - 1)
        A = laplacian_plus_I_operator(D2x, D2y)
        total = size(A, 1)
        rhs = fill(3.0, total)
        mask = falses(ny, nx)
        mask[2:end-1, 2:end-1] .= true
        boundary_values = fill(1.0, ny, nx)
        Ared, bred, boundary_vec = apply_dirichlet!(A, rhs, mask, boundary_values)
        nint = count(identity, mask)
        @test size(Ared, 1) == nint
        @test length(bred) == nint
        @test length(boundary_vec) == total
        sol = solve_poisson_like(Ared, bred)
        full = rebuild_solution(sol, boundary_values, mask)
        @test size(full) == size(boundary_values)
        @test maximum(abs.(full[.!mask] .- 1.0)) < 1e-12
    end

    @testset "Interpolation" begin
        x = cheb_lobatto_nodes(18)
        f = sin.(x)
        xnew = collect(range(-1, 1; length=25))
        fint = barycentric_interp_1d(x, f, xnew)
        @test maximum(abs.(fint .- sin.(xnew))) < 1e-8

        y = cheb_lobatto_nodes(12)
        X, Y = grid2D(x, y)
        F = sin.(X .+ Y)
        xfine = collect(range(-1, 1; length=15))
        yfine = collect(range(-1, 1; length=20))
        Finterp = interp2D_spectral(x, y, F, xfine, yfine)
        exact = sin.(repeat(reshape(xfine, 1, :), length(yfine), 1) .+
                     repeat(reshape(yfine, :, 1), 1, length(xfine)))
        @test maximum(abs.(Finterp .- exact)) < 5e-8
    end
end
