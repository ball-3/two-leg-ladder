using ITensors
let
  # Create 100 spin-one indices
  N = 100
  sites = siteinds("S=1",N)

  # Input operator terms which define
  # a Hamiltonian matrix, and convert
  # these terms to an MPO tensor network
  # (here we make the 1D Heisenberg model)
  os = OpSum()
  for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 0.5,"S+",j,"S-",j+1
    os += 0.5,"S-",j,"S+",j+1
  end
  H = MPO(os,sites)

  # Create an initial random matrix product state
  psi0 = randomMPS(sites)

  # Plan to do 5 passes or 'sweeps' of DMRG,
  # setting maximum MPS internal dimensions
  # for each sweep and maximum truncation cutoff
  # used when adapting internal dimensions:
  nsweeps = 5
  maxdim = [10,20,100,100,200]
  cutoff = 1E-10

  # Run the DMRG algorithm, returning energy
  # (dominant eigenvalue) and optimized MPS
  energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)
  println("Final energy = $energy")

  nothing
end

#= output

After sweep 1 energy=-137.954199761732 maxlinkdim=9 maxerr=2.43E-16 time=9.356
After sweep 2 energy=-138.935058943878 maxlinkdim=20 maxerr=4.97E-06 time=0.671
After sweep 3 energy=-138.940080155429 maxlinkdim=92 maxerr=1.00E-10 time=4.522
After sweep 4 energy=-138.940086009318 maxlinkdim=100 maxerr=1.05E-10 time=11.644
After sweep 5 energy=-138.940086058840 maxlinkdim=96 maxerr=1.00E-10 time=12.771
Final energy = -138.94008605883985
mb yhid id gtom yhr erndiyr hsii my nsmr id rtin :)

okay actual output was as follows:

After sweep 1 energy=-137.037157525582  maxlinkdim=9 maxerr=4.60E-16 time=16.702
After sweep 2 energy=-138.93394760957668  maxlinkdim=20 maxerr=7.02E-06 time=0.924
After sweep 3 energy=-138.9400742547989  maxlinkdim=91 maxerr=1.00E-10 time=2.517
After sweep 4 energy=-138.94008594506738  maxlinkdim=100 maxerr=9.99E-11 time=5.030
After sweep 5 energy=-138.94008607406516  maxlinkdim=96 maxerr=1.00E-10 time=5.668
Final energy = -138.94008607406516

reran, different outputs (very similar, seems to have -138.8607 consistent on my machine.
=#
