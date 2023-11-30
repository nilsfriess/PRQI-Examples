X = 107.5;
h = 0.01;
N = X * 1 / h;
[A,B,M,x,E,Ep] = assemble_schroedinger(X, N, true, true);

create_spectrum_plot(E, Ep, "spectrum.pdf")
create_spectrum_plot(E, Ep, "spectrum-zoom.pdf", 25, 24, 24);

create_table_results(A+B, M, E, x)
