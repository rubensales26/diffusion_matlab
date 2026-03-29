function run_test()
    addpath("../")
    
    fprintf("Testing test_4nodes... ")
    test_4nodes();
    fprintf("OK\n")
    
    fprintf("Testing test_homo_1EG... ")
    test_homo_1EG();
    fprintf("OK\n")
    
    fprintf("Testing test_het_1EG... ")
    test_het_1EG();
    fprintf("OK\n")
end
