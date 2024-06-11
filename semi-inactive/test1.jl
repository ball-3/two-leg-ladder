using ITensors
let
  i = Index(3)     # Index of dimension 3
  @show dim(i)     # = 3
  @show id(i)      # = 0x5d28aa559dd13001 or similar

  ci = copy(i)
  @show ci == i    # = true

  j = Index(5,"j") # Index with a tag "j"

  @show j == i     # = false

  s = Index(2,"n=1,Site") # Index with two tags,
                          # "Site" and "n=1"
  @show hastags(s,"Site") # = true
  @show hastags(s,"n=1")  # = true

  i1 = prime(i) # i1 has a "prime level" of 1
                # but otherwise same properties as i
  @show i1 == i # = false, prime levels do not match

  nothing
end
#=
Test Successful! Output was as follows:

dim(i) = 3
id(i) = 0x5ade15a78b9093bf
ci == i = true
j == i = false
hastags(s, "Site") = true
hastags(s, "n=1") = true
i1 == i = false

=#
