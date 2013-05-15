LIBDIR=/corral-repl/utexas/poldracklab/software_lonestar/local/lib/python

install:
	cp fmriqa.py /corral-repl/utexas/poldracklab/software_lonestar/local/bin/
	cp compute_fd.py $(LIBDIR)
	cp MAD.py $(LIBDIR)
	cp mk_report.py $(LIBDIR)
	cp mk_slice_mosaic.py $(LIBDIR)
	cp plot_timeseries.py $(LIBDIR)