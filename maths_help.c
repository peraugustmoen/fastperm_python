int square(int num) {
    return num * num;
}

void square_vec(double * vec, int n) {
    for (int i = 0; i < n; ++i)
    {
        vec[i] = vec[i]*vec[i];
    }
}