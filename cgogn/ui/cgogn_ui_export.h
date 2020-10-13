
#ifndef CGOGN_UI_EXPORT_H
#define CGOGN_UI_EXPORT_H

#ifdef CGOGN_UI_STATIC_DEFINE
#  define CGOGN_UI_EXPORT
#  define CGOGN_UI_NO_EXPORT
#else
#  ifndef CGOGN_UI_EXPORT
#    ifdef cgogn_ui_EXPORTS
        /* We are building this library */
#      define CGOGN_UI_EXPORT 
#    else
        /* We are using this library */
#      define CGOGN_UI_EXPORT 
#    endif
#  endif

#  ifndef CGOGN_UI_NO_EXPORT
#    define CGOGN_UI_NO_EXPORT 
#  endif
#endif

#ifndef CGOGN_UI_DEPRECATED
#  define CGOGN_UI_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CGOGN_UI_DEPRECATED_EXPORT
#  define CGOGN_UI_DEPRECATED_EXPORT CGOGN_UI_EXPORT CGOGN_UI_DEPRECATED
#endif

#ifndef CGOGN_UI_DEPRECATED_NO_EXPORT
#  define CGOGN_UI_DEPRECATED_NO_EXPORT CGOGN_UI_NO_EXPORT CGOGN_UI_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CGOGN_UI_NO_DEPRECATED
#    define CGOGN_UI_NO_DEPRECATED
#  endif
#endif

#endif
