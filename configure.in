dnl Process this file with autoconf to produce a configure script.
AC_INIT(src/biginteger.cc)


dnl People often put stuff in /usr/local.
if test /usr/local != "$exec_prefix" -a /usr/local != "$prefix" ; then
  CPPFLAGS="$CPPFLAGS -I/usr/local/include"
  LDFLAGS="$LDFLAGS -L/usr/local/lib"
fi


dnl Check for lib gmp
AC_CHECK_LIB(gmp,__gmpz_ui_sub,, [AC_MSG_ERROR([GNU MP not found, or not 4.1.4 or up, see http://gmplib.org])])



#PKG_CPPFLAGS="${PKG_CPPFLAGS} -I/usr/local/include/ "
#PKG_LIBS="${PKG_LIBS} -L/usr/local/lib/ -lgmp"

PKG_CPPFLAGS="${PKG_CPPFLAGS} -I/usr/local/include/ "
PKG_LIBS="${PKG_LIBS} -L/usr/local/lib/ -lgmp"


 


AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_LIBS)
AC_OUTPUT(src/Makevars)




