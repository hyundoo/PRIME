#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::mat prime_core(arma::mat x, arma::mat C, arma::mat sc_data, const int alp){
    // set the parameters
    const int n_genes = sc_data.n_rows;
    const int n_cells = sc_data.n_cols;
    int nbr_len = 0;
    
    // variables
    arma::mat M = x;
    arma::mat sppr_init(n_cells, 1);
    arma::uvec scale_ids;
    arma::uvec ids;
    arma::uvec nbr_id;
    arma::mat sorted_ppr(n_cells, 1);
    arma::mat ppr_app(n_cells, 1);
    arma::mat ppr_app2(n_genes, 1);
    arma::mat ppr_app3(n_genes, 1);
    arma::mat mx(n_genes, 1);   
    arma::mat vx(n_genes, 1);   
    arma::mat scale(n_genes, 1);
    arma::mat g_val(n_genes, 1);
    arma::mat dist_px(n_genes, 1);
    arma::mat imputed_data(n_genes, n_cells);
    
    arma::mat center;
    imputed_data.zeros();
    
    for (int i=0; i<n_cells; i++){
       
        sppr_init.zeros();
        sppr_init(i) = 1;

        ppr_app = (M * sppr_init);
        ppr_app2 = (M * ppr_app);
        ppr_app3 = (M * ppr_app2);
        
        ids = find(C.row(i) >= 0.5);
        nbr_len = ceil(1.25*ids.n_rows);
        nbr_len = min(NumericVector::create(nbr_len, ceil(0.2*n_cells)));

        sorted_ppr = arma::sort(ppr_app3, "descend");
        nbr_id = find(ppr_app3 >=  sorted_ppr(nbr_len));
        
        arma::mat sc_sub = sc_data.cols(nbr_id);
        mx = mean(sc_sub , 1);
        vx = (var(sc_sub , 0, 1) + 0.001);

        g_val = sc_data.col(i);

        scale = 2*exp(-0.05*vx);
        scale_ids = find(scale <= 1);
        scale.rows(scale_ids) = ones(scale_ids.n_rows,1)*1;
                
        dist_px = 1/(1+exp(-(scale %(g_val - mx))));
        imputed_data.col(i) = (dist_px % g_val) + ((1-dist_px) % mx);
    }
    return(imputed_data);
}
