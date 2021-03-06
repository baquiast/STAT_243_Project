library('RUnit')

source('1.R')

test.suite <- defineTestSuite("sub functions",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)