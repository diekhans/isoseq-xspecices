root = ..
include ${root}/defs.mk

build: build_src

build_src: ${DBS:%=%_build_src}
%_build_src:
	${MAKE} build_dest srcDb=$*

build_dest: ${DBS:%=%_build_dest}
%_build_dest: ${srcGp} ${srcMeta} ${srcPsl} ${srcCds} ${srcFa} ${hubSrcBigPsl}
	${MAKE} build_pair srcDb=${srcDb} destDb=$*


ifeq (${srcDb},${destDb})
build_pair:  # does nothing
else
##
# mapping chains
##
build_pair: ${synChains}

${synChains}:
	@mkdir -p $(dir $@)
	buildMappingChains ${srcDb} ${destDb} $@

endif
