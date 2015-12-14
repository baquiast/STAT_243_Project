factorial <- function(n)
{ if (n == 0) { return(1) }
  else { return(n * factorial(n - 1))}
}

test.lower <- function() {
  checkEquals(6, factorial(3))
  checkEqualsNumeric(6, factorial(3))
  checkIdentical(6, factorial(3))
  checkTrue(2 + 2 == 4, 'Arithmetic works')
  checkException(log('a'), 'Unable to take the log() of a string')
}

test.upper <- function() {
  checkEquals(6, factorial(3))
  checkEqualsNumeric(6, factorial(3))
  checkIdentical(6, factorial(3))
  checkTrue(2 + 2 == 4, 'Arithmetic works')
  checkException(log('a'), 'Unable to take the log() of a string')
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}