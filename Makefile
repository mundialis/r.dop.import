MODULE_TOPDIR = ../..

PGM = r.dop.import

ETCFILES = download_urls

include $(MODULE_TOPDIR)/include/Make/Python.make
include $(MODULE_TOPDIR)/include/Make/Script.make

default: script
