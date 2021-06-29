#!/bin/bash

source ./export.sh ; python -c "import database_helpers as dbh;dbh.get_that_there_to_local('$1')"