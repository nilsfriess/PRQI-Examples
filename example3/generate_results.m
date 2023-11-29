[A,B,M,x,E,Ep] = assemble_schroedinger(107.5, 10750, true, true);

create_spectrum_plot(E, Ep, "spectrum.pdf")
create_spectrum_plot(E, Ep, "spectrum-zoom.pdf", 25, 24, 24);

create_table_results(A+B, M, E, x)
