import itertools as it

def get_mean_difference(x, y, debug=False):

    if len(x) != len(y):
        raise ValueError("Please mention the same x, y list")

    total_order = len(x)
    for order in range(1, total_order):
        y = [(y[i]-y[i-1])/(x[i]-x[i-order]) if i>=order else y[i] for i in range(total_order)]
        if debug:
            print("Order:{}\t Caculating...{:>40}{sign}".format(order,repr(y), sign="(FINAL ORDER)" if order==total_order-1 else ""))
    
    return y

class BuildNewtonFunc():
    def __init__(self, inter_points, coff):
        self.inter_points = inter_points
        self.coff = coff
    
    def __str__(self):
        def make_sub_pa(x, variable_name='x'):
            return "({}-{})".format(variable_name, x)

        self.sub_pa_list = [make_sub_pa(x) for x in self.inter_points[:-1]]
        self.str_coff = [str(coff) for coff in self.coff]
        sep_list = ["*".join((self.str_coff[i], *self.sub_pa_list[:i])) for i in range(len(self.str_coff))]
        # print(sep_list)
        
        full_str = ' + '.join(sep_list)

        return full_str
    
    def __call__(self):
        def build_func(x=0):
            p = self.coff[-1]


            #TO-DO: use zip, reversed, etc
            #print(list(zip(self.inter_points, self.coff))[:-1]))
            # for a, x_p in zip(reversed(self.coff[:-1]), reversed(self.coff[:-1])):
            #     p = p*(x-x_p)+a
            for i in range(len(self.inter_points)-1, -1, -1):
                p = p*(x-self.inter_points[i])+self.coff[i]
            return p 
        
        def apply_func(iterobj):
            try:
                for obj in iterobj:
                    yield build_func(obj)
            except StopIteration:
                yield None 
        
        return apply_func


if __name__ == "__main__":
    x = [0, 2, 3, 5]
    y = [1, 3, 2, 5]
    a = get_mean_difference(x, y, debug=True)

    newton_model = BuildNewtonFunc(x, a)
    print("The Newton-interpolation function is:", newton_model)
    func = newton_model()
    var_x = [7]
    for i, j in zip(var_x, func(var_x)):
        print("f({})={}".format(i, j))

    




