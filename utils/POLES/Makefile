
all: lib_brodeau.so drown_kara.so

lib_brodeau.so:
	f2py -c -m lib_brodeau lib_brodeau.f90

drown_kara.so:
	f2py -c -m drown_kara drown_kara.f90

clean:
	@rm lib_brodeau.so drown_kara.so
