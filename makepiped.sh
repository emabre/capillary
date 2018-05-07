#!/bin/bash
make &> make.out
# now you can search for warnings/errors inside the output of make.
# To fin the errors/warning only concerning the files in your directory,
# use the regex: ""[A-Za-z0-9]+\.c" with no quotes, beware the space at the beginning!
