Support for arrays is currently broken.

Fixed-length arrays don't break the compiler, but do have trouble with:

1. Allocating intermediates correctly. Something is wrong with the code that identifies whether or not an intermediate value is an array type.
2. Evaluating special functions on arrays. Things like `log` and `pow` get passed through math_lib, which doesn't support them elementally.

-Adam Jermyn January 25 2021