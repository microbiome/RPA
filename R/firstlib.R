.First.lib <- function(lib, pkg)
{
   library.dynam('RPA', pkg, lib)
   #library.dynam('../src/netresponse', pkg, lib)
   cat('RPA Copyright (C) 2008-2010 Leo Lahti. This program comes with ABSOLUTELY NO WARRANTY.  This is free software, and you are welcome to redistribute it under GNU GPL >=2, see the licensing terms for details.')
    }
