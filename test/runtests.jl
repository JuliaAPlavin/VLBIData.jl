using TestItems
using TestItemRunner
@run_package_tests



@testitem "_" begin
    import Aqua
    Aqua.test_all(VLBIData; ambiguities=true, piracies=true)

    import CompatHelperLocal as CHL
    CHL.@check()
end
