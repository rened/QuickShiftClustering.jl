using FunctionalData, QuickShiftClustering
data = @p map unstack(1:10) (x->10*randn(2,1).+randn(2,1000)) | flatten
@time a = quickshift(data)
Profile.clear_malloc_data()
@time a = quickshift(data)
exit()
