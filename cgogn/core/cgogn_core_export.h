
#ifndef CGOGN_CORE_EXPORT_H
#define CGOGN_CORE_EXPORT_H

#ifdef CGOGN_CORE_STATIC_DEFINE
#  define CGOGN_CORE_EXPORT
#  define CGOGN_CORE_NO_EXPORT
#else
#  ifndef CGOGN_CORE_EXPORT
#    ifdef cgogn_core_EXPORTS
        /* We are building this library */
#      define CGOGN_CORE_EXPORT 
#    else
        /* We are using this library */
#      define CGOGN_CORE_EXPORT 
#    endif
#  endif

#  ifndef CGOGN_CORE_NO_EXPORT
#    define CGOGN_CORE_NO_EXPORT 
#  endif
#endif

#ifndef CGOGN_CORE_DEPRECATED
#  define CGOGN_CORE_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CGOGN_CORE_DEPRECATED_EXPORT
#  define CGOGN_CORE_DEPRECATED_EXPORT CGOGN_CORE_EXPORT CGOGN_CORE_DEPRECATED
#endif

#ifndef CGOGN_CORE_DEPRECATED_NO_EXPORT
#  define CGOGN_CORE_DEPRECATED_NO_EXPORT CGOGN_CORE_NO_EXPORT CGOGN_CORE_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CGOGN_CORE_NO_DEPRECATED
#    define CGOGN_CORE_NO_DEPRECATED
#  endif
#endif

#endif
