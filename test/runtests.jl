using Preferences: set_preferences!
set_preferences!("VLBIData", "dispatch_doctor_mode" => "error")

using TestItems
using TestItemRunner
@run_package_tests

@testitem "Antenna" begin
    using AccessorsExtra

    a = Antenna(:MyAnt, xyz=(1.0, 2.0, 3.0))
    b = Antenna(:XXAnt, xyz=(2, 3, 5))
    @test antenna_names(a) == (:MyAnt,)
    @test antenna_names((a, b)) == (:MyAnt, :XXAnt)
    @test antenna_names([a, b]) == [:MyAnt, :XXAnt]
    @test (@oget UV(a) UV((a, b)) 123) == 123
    
    @test Antenna(:A) < Antenna(:B)
    @test Antenna(:A) == Antenna(:A)

    @test string(a) == "Antenna MyAnt"
end

@testitem "lonlatalt" begin
    using StaticArrays
    A = 6378137.0; f = 1 / 298.257223563; B = A * (1 - f); e2 = f * (2 - f)
    # analytic anchors pin the ellipsoid constants
    g0 = lonlatalt(Antenna(:E, xyz=(A + 100, 0, 0)))        # equator on the prime meridian, 100 m up
    @test g0.lon ≈ 0 atol=1e-12
    @test g0.lat ≈ 0 atol=1e-9
    @test g0.alt ≈ 100 atol=1e-6
    gp = lonlatalt(Antenna(:P, xyz=(0, 0, B + 50)))         # north pole, 50 m up (division-free alt)
    @test gp.lat ≈ π/2 atol=1e-9
    @test gp.alt ≈ 50 atol=1e-6
    # round-trip: geodetic → ECEF (forward) → lonlatalt must invert it
    lon0, lat0, alt0 = deg2rad(-108.12), deg2rad(34.30), 2365.0
    N = A / sqrt(1 - e2 * sin(lat0)^2)
    xyz = SVector((N + alt0) * cos(lat0) * cos(lon0), (N + alt0) * cos(lat0) * sin(lon0), (N * (1 - e2) + alt0) * sin(lat0))
    g = lonlatalt(Antenna(:R, xyz=xyz))
    @test g.lon ≈ lon0 atol=1e-10
    @test g.lat ≈ lat0 atol=1e-10
    @test g.alt ≈ alt0 atol=1e-4
    @test abs(g.lat) < π/2                                   # bare RADIANS, not degrees
end

@testitem "Baseline" begin
    using AccessorsExtra
    using StaticArrays

    bl = Baseline((2, 3))
    @test antenna_names(bl) == (:ANT2, :ANT3)
    @test bl == Baseline((2, 3))
    @test bl != Baseline((:ANT2, :ANT3))
    @test bl != Baseline((2, 4))
    @test bl != Baseline((:ANT2, :ANT4))

    bl = Baseline((:A, :B))
    @test antenna_names(bl) == (:A, :B)
    @test bl == Baseline((:A, :B))
    @test bl != Baseline((:A, :C))

    @test antenna_names(bl) == (:A, :B)
    Accessors.test_getset_laws(antenna_names, bl, (:XX, :YY), (:ZZ, :WW))

    @test (@oget UV(bl) 123) == 123

    @test Baseline((;spec=bl)) == bl
    @test antenna_names((;spec=bl)) == antenna_names(bl)

    @test string(bl) == "Baseline A - B"
end

@testitem "UV" begin
    using AccessorsExtra

    @test (@oget UV([1,2]) 123) == UV(1, 2)

    uv = UV(1, 2)
    @test (@oget UV(uv) 123) == uv
    @test (@oget antenna_names(uv) 123) == 123
    @test (@oget Baseline(uv) 123) == 123

    @test UV((;spec=uv)) == uv
end

@testitem "visspec + visampspec" begin
    using AccessorsExtra
    using StaticArrays

    spec = VLBIData.VisSpec0(UV(10, -20))
    @test UV(spec) == UV(10, -20)
    @test !hasoptic(spec, Baseline)
    @test hasoptic(spec, UV)
    @test !hasoptic(spec, antenna_names)
    @test_broken conj(spec) == VLBIData.VisSpec0(UV(-10, 20))
    @test string(spec) == "Vis, 22.4"
    @test spec == VLBIData.VisSpec0(UV(10, -20.))
    @test spec != VLBIData.VisSpec0(UV(10, -21))

    spec = VisSpec(Baseline((2, 3)), UV(10, -20))
    @test Baseline(spec) == Baseline((2, 3))
    @test UV(spec) == UV(10, -20)
    @test antenna_names(spec) == (:ANT2, :ANT3)
    @test hasoptic(spec, Baseline)
    @test hasoptic(spec, UV)
    @test hasoptic(spec, antenna_names)
    @test conj(spec) == VisSpec(Baseline((3, 2)), UV(-10, 20))
    @test string(spec) == "Vis ANT2 - ANT3, 22.4"
    @test spec == VisSpec(Baseline((2, 3)), UV(10, -20.))
    @test spec != VisSpec(Baseline((2, 4)), UV(10, -20))
    @test spec != VisSpec(Baseline((2, 3)), UV(10, -21))

    aspec = VisAmpSpec(spec)
    @test aspec == VisAmpSpec(Baseline((2, 3)), UV(10, -20))
    @test Baseline(aspec) == Baseline((2, 3))
    @test UV(aspec) == UV(10, -20)
    @test antenna_names(aspec) == (:ANT2, :ANT3)
    @test VisSpec(aspec) == spec
    @test conj(aspec) == VisAmpSpec(VisSpec(Baseline((3, 2)), UV(-10, 20)))
    @test string(aspec) == "VisAmp ANT2 - ANT3, 22.4"
    @test aspec == VisAmpSpec(Baseline((2, 3)), UV(10, -20.))
    @test aspec != VisAmpSpec(Baseline((2, 4)), UV(10, -20))
    @test aspec != VisAmpSpec(Baseline((2, 3)), UV(10, -21))

    @test VLBIData.AbstractSpec((;spec)) == spec
    @test VisSpec((;spec)) == spec
    @test_throws "expected V" VisAmpSpec((;spec))
    @test Baseline((;spec)) == Baseline(spec)
    @test antenna_names((;spec)) == antenna_names(spec)
    @test UV((;spec)) == UV(spec)
    @test VLBIData.conjvis((;spec, value=1+2im)) == (spec=conj(spec), value=1-2im)
    @test VLBIData.conjvis((;spec, value=1+2im, stokes=:RR)) == (spec=conj(spec), value=1-2im, stokes=:RR)
    @test VLBIData.conjvis((;spec, value=1+2im, stokes=:RL)) == (spec=conj(spec), value=1-2im, stokes=:LR)
    @test VLBIData.conjvis((;spec)) == (spec=conj(spec),)
    @test VLBIData.conjvis((;spec, stokes=:RL)) == (spec=conj(spec), stokes=:LR)
end

@testitem "frequencies" begin
    using Unitful
    struct MyFreq
        freq::Float64
    end
    VLBIData.frequency(x::MyFreq) = x.freq
    @test frequency((;freq_spec=MyFreq(10))) == 10
    @test frequency(10u"GHz") == 10u"GHz"
    @test frequency((;freq_spec=10u"GHz")) == 10u"GHz"
    @test_throws MethodError frequency(10u"m")
end

@testitem "renumbering" begin
    uvtbl_orig = [
        (spec=VLBIData.VisSpec0(UV(0, 1)), freq_spec=123, stokes=:RR, value=1),
        (spec=VLBIData.VisSpec0(UV(0, 2)), freq_spec=123, stokes=:RR, value=2),
        (spec=VLBIData.VisSpec0(UV(0, 5)), freq_spec=124, stokes=:LL, value=3),
    ]
    @test VLBIData.uv_reindex(uvtbl_orig) == [
        (spec=VLBIData.VisSpec0(UV(0, 1)), freq_spec=123, stokes=:RR, value=1),
        (spec=VLBIData.VisSpec0(UV(0, 2)), freq_spec=123, stokes=:RR, value=2),
        (spec=VLBIData.VisSpec0(UV(0, 5)), freq_spec=124, stokes=:LL, value=3),
    ]

    uvtbl_orig = [
        (spec=VisSpec(Baseline((:A, :B)), UV(0, 1)), stokes=:RR, value=1),
        (spec=VisSpec(Baseline((:A, :B)), UV(0, 2)), stokes=:RR, value=2),
        (spec=VisSpec(Baseline((:C, :B)), UV(0, 5)), stokes=:LL, value=3),
    ]
    @test VLBIData.uv_reindex(uvtbl_orig) == uvtbl_orig

    uvtbl_orig = [
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 1))), freq_spec=123, stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 2))), freq_spec=123, stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline((:C, :B)), UV(0, 5))), freq_spec=124, stokes=:LL, value=3),
    ]
    @test VLBIData.uv_reindex(uvtbl_orig) == uvtbl_orig

    uvtbl_orig = [
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 1))), freq_spec=(;ix=123), stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 2))), freq_spec=(;ix=123), stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline((:C, :B)), UV(0, 5))), freq_spec=(;ix=124), stokes=:LL, value=3),
    ]
    @test VLBIData.uv_reindex(uvtbl_orig) == [
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 1))), freq_spec=(;ix=1), stokes=:RR, value=1),
        (spec=VisAmpSpec(VisSpec(Baseline((:A, :B)), UV(0, 2))), freq_spec=(;ix=1), stokes=:RR, value=2),
        (spec=VisAmpSpec(VisSpec(Baseline((:C, :B)), UV(0, 5))), freq_spec=(;ix=2), stokes=:LL, value=3),
    ]
end

@testitem "comradebase" begin
    import ComradeBase as CB
    using Unitful

    p = (; U = 0.1 * randn(60), V = 0.1 * randn(60), Ti = collect(Float64, 1:60), Fr = fill(230.0e9, 60))
    d = rand(60)
    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p))
    uvtbl = uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :freq_spec, :value, :datetime))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].freq_spec ≈ 230u"GHz"
    @test uvtbl[10].value == d[10]
    @test uvtbl[10].datetime == p.Ti[10]

    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p[(:U, :V)]))
    uvtbl = uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :value))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].value == d[10]

    m = CB.UnstructuredMap(d, CB.UnstructuredDomain(p[(:U, :V, :Fr)]))
    uvtbl = uvtable(m)
    @test length(uvtbl) == 60
    @test issetequal(keys(uvtbl[10]), (:spec, :value, :freq_spec))
    @test uvtbl[10].spec::UV == UV(p.U[10], p.V[10])
    @test uvtbl[10].value == d[10]
    @test uvtbl[10].freq_spec ≈ 230u"GHz"
end

@testitem "_parallactic_angle" begin
    using StaticArrays

    # Closed-form standard parallactic angle (geocentric latitude from xyz):
    # PA = atan2(sin(H_local)*cos(lat), cos(dec)*sin(lat) - sin(dec)*cos(lat)*cos(H_local))
    # where H_local = H + lon, lat = atan2(z, sqrt(x²+y²)), lon = atan2(y, x)
    function par_standard(xyz, H, dec)
        x, y, z = xyz
        r_xy = hypot(x, y)
        sin_lat, cos_lat = z / hypot(r_xy, z), r_xy / hypot(r_xy, z)
        H_local = H + atan(y, x)
        atan(sin(H_local)*cos_lat, cos(dec)*sin_lat - sin(dec)*cos_lat*cos(H_local))
    end

    cases = [
        (SVector(6.378e6 * cos(deg2rad(45)), 0.0, 6.378e6 * sin(deg2rad(45))), 0.0, 0.0),
        (SVector(6.378e6, 0.0, 0.0), -π/2, 0.0),
        (SVector(6.378e6 * cos(deg2rad(45)) * cos(deg2rad(30)),
                 6.378e6 * cos(deg2rad(45)) * sin(deg2rad(30)),
                 6.378e6 * sin(deg2rad(45))), deg2rad(15), deg2rad(20)),
        (SVector(6.378e6 * cos(deg2rad(-30)) * cos(deg2rad(120)),
                 6.378e6 * cos(deg2rad(-30)) * sin(deg2rad(120)),
                 6.378e6 * sin(deg2rad(-30))), deg2rad(-90), deg2rad(-45)),
    ]
    for (xyz, H, dec) in cases
        @test VLBIData._parallactic_angle(xyz, H, dec) ≈ par_standard(xyz, H, dec) atol=1e-12
    end

    # Hardcoded reference: VLBA-style antenna at lat=45°N, lon=0°, source at meridian (H=0), dec=0 → PA=0
    xyz_eq_meridian = SVector(6.378e6 * cos(deg2rad(45)), 0.0, 6.378e6 * sin(deg2rad(45)))
    @test VLBIData._parallactic_angle(xyz_eq_meridian, 0.0, 0.0) ≈ 0.0 atol=1e-12

    # Hardcoded reference: equator antenna, source rising in east (H=-π/2), dec=0 → PA=-π/2
    xyz_eq = SVector(6.378e6, 0.0, 0.0)
    @test VLBIData._parallactic_angle(xyz_eq, -π/2, 0.0) ≈ -π/2 atol=1e-12

    # Hardcoded reference value at a non-trivial configuration:
    # lat=45°N, lon=30°E, dec=20°, H=15° → PA from independent calc
    xyz_nontriv = SVector(6.378e6 * cos(deg2rad(45)) * cos(deg2rad(30)),
                          6.378e6 * cos(deg2rad(45)) * sin(deg2rad(30)),
                          6.378e6 * sin(deg2rad(45)))
    @test rad2deg(VLBIData._parallactic_angle(xyz_nontriv, deg2rad(15), deg2rad(20))) ≈ 45.377584812297336 atol=1e-10
end

@testitem "polarizedtypes extension (manual table)" begin
    using VLBIData

    spec1 = VisSpec(Baseline((:A, :B)), UV(1.0, 2.0))
    spec2 = VisSpec(Baseline((:C, :D)), UV(3.0, 4.0))

    # Create a minimal fake uvtable with all 4 stokes
    uvtbl = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RL, value=13.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=14.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:RR, value=21.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:LR, value=22.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:RL, value=23.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:LL, value=24.0),
    ]
    cm_tbl = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl)

    @test length(cm_tbl) == 2
    @test cm_tbl[1].value isa CoherencyMatrix
    @test cm_tbl[2].value isa CoherencyMatrix
    @test cm_tbl[1].value == [11 13; 12 14]
    @test cm_tbl[2].value == [21 23; 22 24]
    @test !haskey(cm_tbl[1], :stokes)

    # Partial group (missing :LL) → CoherencyMatrix with NaN in the LL slot
    uvtbl_missing_LL = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RL, value=13.0),
    ]
    cm_missing_LL = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_missing_LL)
    @test length(cm_missing_LL) == 1
    @test cm_missing_LL[1].value[1, 1] == 11.0  # :RR
    @test cm_missing_LL[1].value[2, 1] == 12.0  # :LR
    @test cm_missing_LL[1].value[1, 2] == 13.0  # :RL
    @test isnan(cm_missing_LL[1].value[2, 2])   # :LL filled with NaN

    # IPol: from all 4 stokes (picks parallel hands)
    ipol_tbl = VLBIData.uvtable_values_to(IPol, uvtbl)
    @test length(ipol_tbl) == 2
    @test ipol_tbl[1].value == (11.0 + 14.0) / 2
    @test ipol_tbl[2].value == (21.0 + 24.0) / 2
    @test !haskey(ipol_tbl[1], :stokes)

    # IPol: from only parallel hands
    uvtbl_par = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=10.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=20.0),
    ]
    ipol_par = VLBIData.uvtable_values_to(IPol, uvtbl_par)
    @test length(ipol_par) == 1
    @test ipol_par[1].value == 15.0

    # IPol: linear basis
    uvtbl_lin = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:XX, value=6.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:YY, value=4.0),
    ]
    ipol_lin = VLBIData.uvtable_values_to(IPol, uvtbl_lin)
    @test ipol_lin[1].value == 5.0

    # IPol: from single parallel hand (RR only)
    uvtbl_single_rr = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RL, value=13.0),
    ]
    ipol_single = VLBIData.uvtable_values_to(IPol, uvtbl_single_rr)
    @test length(ipol_single) == 1
    @test ipol_single[1].value == 11.0
    @test !haskey(ipol_single[1], :stokes)

    # IPol: from single parallel hand (LL only)
    uvtbl_single_ll = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=7.0),
    ]
    ipol_single_ll = VLBIData.uvtable_values_to(IPol, uvtbl_single_ll)
    @test length(ipol_single_ll) == 1
    @test ipol_single_ll[1].value == 7.0

    # IPol: when no parallel hands
    uvtbl_no_par = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RL, value=13.0),
    ]
    @test isempty(VLBIData.uvtable_values_to(IPol, uvtbl_no_par))

    # IPol: data already has stokes=:I (passthrough)
    uvtbl_stokesI = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:I, value=42.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:I, value=99.0),
    ]
    ipol_from_I = VLBIData.uvtable_values_to(IPol, uvtbl_stokesI)
    @test length(ipol_from_I) == 2
    @test ipol_from_I[1].value == 42.0
    @test ipol_from_I[2].value == 99.0
    @test !haskey(ipol_from_I[1], :stokes)

    # Same baseline, different UV per stokes: should still group together
    bl1 = Baseline((:A, :B))
    bl2 = Baseline((:C, :D))
    uvtbl_diffuv = [
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.0, 2.0)),  stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.01, 2.01)), stokes=:LR, value=12.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.02, 2.02)), stokes=:RL, value=13.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.03, 2.03)), stokes=:LL, value=14.0),
    ]
    cm_diffuv = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_diffuv)
    @test length(cm_diffuv) == 1
    @test cm_diffuv[1].value == [11 13; 12 14]

    ipol_diffuv = VLBIData.uvtable_values_to(IPol, uvtbl_diffuv)
    @test length(ipol_diffuv) == 1
    @test ipol_diffuv[1].value == (11.0 + 14.0) / 2
    @test UV(ipol_diffuv[1]) ≈ [1.015, 2.015]

    # Multiple baselines with differing UV: groups should be per-baseline
    uvtbl_multi = [
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.0, 2.0)),  stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl1, UV(1.01, 2.01)), stokes=:LL, value=14.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl2, UV(3.0, 4.0)),  stokes=:RR, value=31.0),
        (datetime=1, freq_spec=100.0, spec=VisSpec(bl2, UV(3.01, 4.01)), stokes=:LL, value=34.0),
    ]
    ipol_multi = VLBIData.uvtable_values_to(IPol, uvtbl_multi)
    @test length(ipol_multi) == 2
    vals = sort([r.value for r in ipol_multi])
    @test vals == [(11.0 + 14.0) / 2, (31.0 + 34.0) / 2]

    # CoherencyMatrix: partial stokes per group → fill missing slots with NaN
    uvtbl_partial = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=14.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:RR, value=21.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:LR, value=22.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:RL, value=23.0),
        (datetime=2, freq_spec=200.0, spec=spec2, stokes=:LL, value=24.0),
    ]
    cm_partial = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_partial)
    @test length(cm_partial) == 2
    # bl1: only :RR (slot [1,1]) and :LL (slot [2,2]) present; cross-hands are NaN
    @test cm_partial[1].value[1, 1] == 11.0
    @test cm_partial[1].value[2, 2] == 14.0
    @test isnan(cm_partial[1].value[2, 1])  # :LR slot
    @test isnan(cm_partial[1].value[1, 2])  # :RL slot
    # bl2: complete — all four finite
    @test cm_partial[2].value == [21 23; 22 24]

    # ComplexF32 input must stay ComplexF32 (no Float64 promotion via NaN fallback)
    uvtbl_partial_f32 = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=ComplexF32(11.0f0)),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=ComplexF32(14.0f0)),
    ]
    cm_partial_f32 = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_partial_f32)
    @test eltype(cm_partial_f32[1].value) === ComplexF32

    # Parallel-hands-only group (:RR, :LL) → cross-hands NaN, parallel hands finite
    uvtbl_par_only = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=10.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=20.0),
    ]
    cm_par_only = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_par_only)
    @test length(cm_par_only) == 1
    @test cm_par_only[1].value[1, 1] == 10.0
    @test cm_par_only[1].value[2, 2] == 20.0
    @test isnan(cm_par_only[1].value[2, 1])
    @test isnan(cm_par_only[1].value[1, 2])
    @test !haskey(cm_par_only[1], :stokes)

    # Duplicate stokes within a group: the affected slot becomes NaN
    # (hasoptic for filteronly is `count == 1`, so 2+ matches → fallback)
    uvtbl_dup_RR = [
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=11.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:RR, value=99.0),
        (datetime=1, freq_spec=100.0, spec=spec1, stokes=:LL, value=14.0),
    ]
    cm_dup = VLBIData.uvtable_values_to(CoherencyMatrix, uvtbl_dup_RR)
    @test length(cm_dup) == 1
    @test isnan(cm_dup[1].value[1, 1])  # :RR duplicated → NaN
    @test cm_dup[1].value[2, 2] == 14.0 # :LL still finite
end

@testitem "IPol combined row: uv, count" begin
    using VLBIData
    using AccessorsExtra
    using Uncertain
    using Dates

    bl = Baseline((:A, :B))
    t = DateTime(2020, 1, 1, 12, 34, 56, 789)
    rr = (datetime=t, freq_spec=100.0, spec=VisSpec(bl, UV(1000.0, -2000.0)), stokes=:RR, value=1.0 ±ᵤ 0.1, count=3)
    ll = (datetime=t, freq_spec=100.0, spec=VisSpec(bl, UV(1010.0, -1990.0)), stokes=:LL, value=3.0 ±ᵤ 0.2, count=5)

    # two hands with different uv → equal-weight mean of values and uv, counts summed
    for tbl in ([rr, ll], [ll, rr])
        r = only(@inferred VLBIData.uvtable_values_to(IPol, tbl))
        @test keys(r) == (:datetime, :freq_spec, :spec, :value, :count)
        @test r.datetime == t
        @test r.freq_spec == 100.0
        @test r.spec == VisSpec(bl, UV(1005.0, -1995.0))
        @test Uncertain.value(r.value) == 2.0
        @test Uncertain.uncertainty(r.value) ≈ 0.11180339887498948
        @test r.count == 8
    end

    # without a count field
    r = only(VLBIData.uvtable_values_to(IPol, [Base.structdiff(x, NamedTuple{(:count,)}) for x in (rr, ll)]))
    @test keys(r) == (:datetime, :freq_spec, :spec, :value)
    @test r.spec == VisSpec(bl, UV(1005.0, -1995.0))

    # one hand → the row unchanged, without stokes
    @test VLBIData.uvtable_values_to(IPol, [ll]) == [Base.structdiff(ll, NamedTuple{(:stokes,)})]

    # identical uv → uv exactly unchanged
    uv = UV(0.1, 0.7)
    r = only(VLBIData.uvtable_values_to(IPol, [(@set rr.spec = VisSpec(bl, uv)), (@set ll.spec = VisSpec(bl, uv))]))
    @test UV(r) === uv
    @test r.datetime === t
end

@testitem "likelihoods" begin
    using InterferometricModels
    using AccessorsExtra
    using StaticArrays
    using Uncertain

    model = construct(CircularGaussian, flux=>1.0, fwhm_average=>0.1, coords=>SVector(0., 0.1))

    visspecs = (
        VisSpec(Baseline((:A, :B)), UV(1, 2)),
        VisSpec(Baseline((:B, :C)), UV(3, 4)),
        VisSpec(Baseline((:C, :A)), UV(5, 6)),
        VisSpec(Baseline((:A, :D)), UV(7, 8)),
    )

    # --- VisSpec ---
    modvis = visibility(model, visspecs[1])
    x_vis = (spec=visspecs[1], value=modvis ±ᵤ 0.1)
    μ, σ, val = VLBIData.loglike_parts(model, x_vis)
    @test μ ≈ modvis
    @test val ≈ modvis  # data == model
    @test σ ≈ 0.1
    @test isfinite(VLBIData.loglike(model, x_vis))

    # perturbed data has lower loglike
    x_vis_pert = (spec=visspecs[1], value=(modvis + 0.5 + 0.5im) ±ᵤ 0.1)
    @test VLBIData.loglike(model, x_vis) > VLBIData.loglike(model, x_vis_pert)

    # --- VisAmpSpec ---
    aspec = VisAmpSpec(visspecs[1])
    x_amp = (spec=aspec, value=abs(modvis) ±ᵤ 0.05)
    μ_a, σ_a, val_a = VLBIData.loglike_parts(model, x_amp)
    @test μ_a ≈ abs(modvis)
    @test val_a ≈ abs(modvis)
    @test isfinite(VLBIData.loglike(model, x_amp))

    x_amp_pert = (spec=aspec, value=(abs(modvis) + 0.5) ±ᵤ 0.05)
    @test VLBIData.loglike(model, x_amp) > VLBIData.loglike(model, x_amp_pert)

    # high SNR: besselix argument val*μ/σ² can be huge
    m_bright = Point(100.0, SVector(0.0, 0.0))
    x_amp_hsnr = (spec=aspec, value=100.0 ±ᵤ 0.001)
    @test isfinite(VLBIData.loglike(m_bright, x_amp_hsnr))

    # --- ClosurePhaseSpec ---
    cpspec = ClosurePhaseSpec(visspecs[1:3])
    mod_bispec = visibility(model, cpspec)
    x_cp = (spec=cpspec, value=mod_bispec ±ᵤ 0.01)
    μ_cp, σ_cp, val_cp = VLBIData.loglike_parts(model, x_cp)
    @test μ_cp ≈ angle(mod_bispec)
    @test isfinite(VLBIData.loglike(model, x_cp))

    # perturbed closure phase
    x_cp_pert = (spec=cpspec, value=(mod_bispec * cis(0.5)) ±ᵤ 0.01)
    @test VLBIData.loglike(model, x_cp) > VLBIData.loglike(model, x_cp_pert)

    # wrapping: phase near π should work
    x_cp_wrap = (spec=cpspec, value=(-abs(mod_bispec) + 0.0im) ±ᵤ 0.1)
    @test isfinite(VLBIData.loglike(model, x_cp_wrap))

    # --- ClosureAmpSpec ---
    caspec = ClosureAmpSpec(visspecs[1:4])
    mod_camp = visibility(model, caspec)
    x_ca = (spec=caspec, value=abs(mod_camp) ±ᵤ 0.05)
    μ_ca, σ_ca, val_ca = VLBIData.loglike_parts(model, x_ca)
    @test μ_ca ≈ log(abs(mod_camp))
    @test isfinite(VLBIData.loglike(model, x_ca))

    x_ca_pert = (spec=caspec, value=(abs(mod_camp) * 3) ±ᵤ 0.05)
    @test VLBIData.loglike(model, x_ca) > VLBIData.loglike(model, x_ca_pert)

    # --- Cross-check against Distributions.jl ---
    using Distributions: logpdf, Normal, Rician, VonMises

    # VisSpec: two independent Normals on Re and Im
    for x_v in [x_vis, x_vis_pert]
        μ_v, σ_v, val_v = VLBIData.loglike_parts(model, x_v)
        expected = logpdf(Normal(real(μ_v), σ_v), real(val_v)) + logpdf(Normal(imag(μ_v), σ_v), imag(val_v))
        @test VLBIData.loglike(model, x_v) ≈ expected
    end

    # VisAmpSpec: Rician
    # XXX: Distributions.jl Rician delegates to Rmath's dnchisq which is inaccurate for large noncentrality parameters (ncp > ~30).
    # R docs warn: "The code for non-zero ncp is principally intended to be used for moderate values of ncp".
    # Our besselix formula matches BigFloat series and scipy to machine precision.
    let
        μ_a, σ_a, val_a = VLBIData.loglike_parts(model, x_amp)
        @test VLBIData.loglike(model, x_amp) ≈ logpdf(Rician(μ_a, σ_a), val_a)
    end
    let
        μ_a, σ_a, val_a = VLBIData.loglike_parts(model, x_amp_pert)
        @test_broken VLBIData.loglike(model, x_amp_pert) ≈ logpdf(Rician(μ_a, σ_a), val_a)
    end

    # ClosurePhaseSpec: VonMises
    for x_c in [x_cp, x_cp_pert, x_cp_wrap]
        μ_c, σ_c, val_c = VLBIData.loglike_parts(model, x_c)
        expected = logpdf(VonMises(μ_c, inv(σ_c^2)), val_c)
        @test VLBIData.loglike(model, x_c) ≈ expected
    end

    # ClosureAmpSpec: Normal
    for x_l in [x_ca, x_ca_pert]
        μ_l, σ_l, val_l = VLBIData.loglike_parts(model, x_l)
        expected = logpdf(Normal(μ_l, σ_l), val_l)
        @test VLBIData.loglike(model, x_l) ≈ expected
    end
end

@testitem "loglike ClosureAmpSpec zero allocations" begin
    using InterferometricModels, Uncertain, StaticArrays, Unitful, UnitfulAngles

    function _test_loglike_allocs()
        sky = EllipticGaussian(flux=1.0, σ_major=0.5u"mas", ratio_minor_major=0.9,
            pa_major=0.0, coords=(0.0, 0.0) .* u"mas")
        vm = Ref(visibility(sky))
        spec = ClosureAmpSpec((
            VisSpec(Baseline((:A, :B)), UV(SVector(1f6, 2f6))),
            VisSpec(Baseline((:B, :C)), UV(SVector(2f6, 1f6))),
            VisSpec(Baseline((:C, :D)), UV(SVector(1.5f6, 1.5f6))),
            VisSpec(Baseline((:A, :D)), UV(SVector(0.5f6, 2.5f6))),
        ))
        x = Ref((; spec, value=Uncertain.Value(0.8 + 0.0im, 0.1)))
        VLBIData.loglike(vm[], x[])  # compile
        @allocated VLBIData.loglike(vm[], x[])
    end

    @test _test_loglike_allocs() == 0
end

@testitem "coherency linear feeds" begin
    using VLBIData
    using VLBIData: PolT
    LinBasis = PolT.LinBasis
    bl = Baseline((:A,:B))
    rows = [(datetime=1, freq_spec=1, spec=VisSpec(bl, UV(0.0,0.0)), stokes=s, value=v)
            for (s,v) in ((:XX,1.0+0im),(:YY,2.0+0im),(:XY,0.0+0im),(:YX,0.0+0im))]
    out = only(VLBIData.uvtable_values_to(CoherencyMatrix, rows)).value
    @test out isa CoherencyMatrix{<:LinBasis,<:LinBasis}
    @test out[1,1] == 1.0+0im && out[2,2] == 2.0+0im
end

@testitem "feed_slot_map" begin
    using StaticArrays
    using VLBIData: feed_slot_map, feed_basis
    @test feed_slot_map([:RR,:LL,:RL,:LR]) == SMatrix{2,2,Int8}(1, 4, 3, 2)   # [1,1]=RR=1 [2,1]=LR=4 [1,2]=RL=3 [2,2]=LL=2
    @test feed_slot_map([:RR,:LL])          == SMatrix{2,2,Int8}(1, 0, 0, 2)   # cross hands absent → 0
    @test feed_slot_map([:RR])              == SMatrix{2,2,Int8}(1, 0, 0, 0)   # single-pol: basis inferred R→(R,L)
    @test feed_slot_map([:XX])              == SMatrix{2,2,Int8}(1, 0, 0, 0)
    @test feed_slot_map([:XX,:YY,:XY,:YX])  == SMatrix{2,2,Int8}(1, 4, 3, 2)
    @test feed_slot_map([:LR,:RL,:LL,:RR])  == SMatrix{2,2,Int8}(4, 1, 2, 3)   # index looked up in the given order
    @test_throws Exception feed_slot_map([:RR,:XX])     # mixed bases
    @test_throws Exception feed_slot_map([:I,:Q,:U,:V]) # Stokes parameters are not a coherency
    @test feed_basis([:RR]) isa typeof(feed_basis([:RR,:LL,:RL,:LR]))   # both circular
    @test (@inferred feed_slot_map([:RR,:LL,:RL,:LR])) == SMatrix{2,2,Int8}(1, 4, 3, 2)   # type-stable
end

@testitem "GapBasedScans" begin
    using Dates, Unitful, StructArrays
    using VLBIData: GapBasedScans, scan_ids, scan_intervals

    dts = DateTime(2020, 1, 1) .+ Second.([0, 2, 4, 6, 8, 10, 200, 202])
    tbl = StructArray((datetime = dts, source = [:A, :A, :A, :B, :B, :B, :B, :B]))
    ids = [1, 1, 1, 2, 2, 2, 3, 3]
    # back-to-back sources are separate scans even with no gap between them
    @test (@inferred scan_ids(GapBasedScans(), tbl)) == ids
    # one source spanning a sub-min_gap boundary stays one scan
    @test scan_ids(GapBasedScans(min_gap=5u"minute"), tbl) == [1, 1, 1, 2, 2, 2, 2, 2]
    # without a source column only gaps split
    @test scan_ids(GapBasedScans(), StructArray((datetime = dts,))) == [1, 1, 1, 1, 1, 1, 2, 2]
    # ids follow time order and land back on their own rows whatever the row order
    perm = [4, 1, 7, 2, 8, 5, 3, 6]
    @test scan_ids(GapBasedScans(), tbl[perm]) == ids[perm]
    @test length(scan_intervals(GapBasedScans(), tbl)) == 3
    # nodding between two sources
    nod = StructArray((datetime = DateTime(2020, 1, 1) .+ Second.(0:2:8), source = [:A, :A, :B, :B, :A]))
    @test scan_ids(GapBasedScans(), nod) == [1, 1, 2, 2, 3]
    # single row, and a single row per source
    @test scan_ids(GapBasedScans(), tbl[1:1]) == [1]
    @test scan_ids(GapBasedScans(), tbl[[1, 4]]) == [1, 2]
    # sources sharing a timestamp: a table that names no antennas has one stream, and one stream
    # cannot be on two sources at once — refused rather than cut by row order
    tied = StructArray((datetime = fill(DateTime(2020, 1, 1), 4), source = [:B, :A, :B, :A]))
    @test_throws "at one instant 2020-01-01T00:00:00; a subnetting schedule needs antenna columns" scan_ids(GapBasedScans(), tied)
end

@testitem "IntegrationGridScans" begin
    using Dates, Unitful, StructArrays
    using VLBIData: GapBasedScans, IntegrationGridScans, scan_ids, grid_position, on_grid

    t0 = DateTime(2020, 1, 1)
    @test grid_position(t0 + Second(9), t0, 2.0u"s") === 4.5
    @test grid_position(t0 + Millisecond(2010), t0, 2.0f0u"s") ≈ 1.005
    @test grid_position(t0 + Millisecond(81131), t0, 0.131072f0u"s") isa Float64
    @test grid_position(t0 + Millisecond(81131), t0, 0.131072f0u"s") ≈ 81.131 / 0.131072  rtol=1e-6
    @test on_grid(1.005)
    @test !on_grid(1.015)
    @test on_grid(-0.01)

    # three antennas, all baselines at each instant, 2 s integrations
    rows(ts, dt) = StructArray([(datetime = t0 + Millisecond(round(Int, 1000t)), baseline = Baseline(bl), source = :S, int_time = dt)
                                for t in ts for bl in ((:A, :B), (:A, :C), (:B, :C))])
    # the grid jumps by half a step: one scan by gaps, two on the grid
    jump = rows([0, 2, 4, 6, 9, 11, 13], 2.0u"s")
    @test scan_ids(GapBasedScans(), jump) == fill(1, 21)
    @test scan_ids(IntegrationGridScans(), jump) == [fill(1, 12); fill(2, 9)]
    # within the tolerance of the grid: one scan
    @test scan_ids(IntegrationGridScans(), rows([0, 2, 4, 6.01, 8, 10.02], 2.0u"s")) == fill(1, 18)
    # off it: a new run anchored at that row, so the next row is off the new grid too
    @test scan_ids(IntegrationGridScans(), rows([0, 2, 4, 6.03, 8], 2.0u"s")) == [fill(1, 9); fill(2, 3); fill(3, 3)]
    # a different int_time starts a new run, on grid or not
    change = [rows([0, 2, 4, 6], 2.0u"s"); rows([8, 9, 10], 1.0u"s")]
    @test scan_ids(GapBasedScans(), change) == fill(1, 21)
    @test scan_ids(IntegrationGridScans(), change) == [fill(1, 12); fill(2, 9)]
    # gaps and sources still cut
    @test scan_ids(IntegrationGridScans(), rows([0, 2, 100, 102], 2.0u"s")) == [fill(1, 6); fill(2, 6)]
    @test (@inferred scan_ids(IntegrationGridScans(), jump)) == scan_ids(IntegrationGridScans(), jump[end:-1:1])[end:-1:1]
    @test_throws "IntegrationGridScans: antenna A observes S and T" scan_ids(IntegrationGridScans(),
        StructArray([(datetime = t0, baseline = Baseline((:A, :B)), source = :S, int_time = 2.0u"s"),
                     (datetime = t0, baseline = Baseline((:A, :C)), source = :T, int_time = 2.0u"s")]))
end

@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIData; piracies=(;broken=true), ambiguities=(;broken=true))

    import CompatHelperLocal as CHL
    CHL.@check()
end
