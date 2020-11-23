
#ifndef CGOGN_IO_EXPORT_H
#define CGOGN_IO_EXPORT_H

#ifdef CGOGN_IO_STATIC_DEFINE
#  define CGOGN_IO_EXPORT
#  define CGOGN_IO_NO_EXPORT
#else
#  ifndef CGOGN_IO_EXPORT
#    ifdef cgogn_io_EXPORTS
        /* We are building this library */
#      define CGOGN_IO_EXPORT 
#    else
        /* We are using this library */
#      define CGOGN_IO_EXPORT 
#    endif
#  endif

#  ifndef CGOGN_IO_NO_EXPORT
#    define CGOGN_IO_NO_EXPORT 
#  endif
#endif

#ifndef CGOGN_IO_DEPRECATED
#  define CGOGN_IO_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CGOGN_IO_DEPRECATED_EXPORT
#  define CGOGN_IO_DEPRECATED_EXPORT CGOGN_IO_EXPORT CGOGN_IO_DEPRECATED
#endif

#ifndef CGOGN_IO_DEPRECATED_NO_EXPORT
#  define CGOGN_IO_DEPRECATED_NO_EXPORT CGOGN_IO_NO_EXPORT CGOGN_IO_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CGOGN_IO_NO_DEPRECATED
#    define CGOGN_IO_NO_DEPRECATED
#  endif
#endif

#endif /* CGOGN_IO_EXPORT_H */
