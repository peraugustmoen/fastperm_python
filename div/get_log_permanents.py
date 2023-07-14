def get_log_perms(X,a,b):
	# do various checks here
	print(X.shape)
	T = X.shape[0]
	n = X.shape[1]

	#X = np.copy(X.transpose(), order='F')
	X = np.transpose(X)
	X = np.array(X,dtype=np.float64,order='F')
	a = np.array(a,dtype=np.float64)
	b = np.array(b,dtype=np.float64)
	T = int(T)
	n = int(n)

	return fastperm.get_log_permanents(X,a,b,n,T,0)