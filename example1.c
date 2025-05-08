#include <stdio.h>

// Define the maximum sizes for the array (you can adjust this based on your needs)
#define MAX_ITER 5    // Number of iterations (mu + 1)
#define MAX_DEG 4     // Maximum degree for x (e.g., x^3)
#define MAX_ALPHA 8   // Maximum degree for alpha (you may need to adjust this)

// Define the sigma array (initialize with random powers of alpha for demonstration)
uint32_t sigma[MAX_ITER][MAX_DEG][MAX_ALPHA];

// Example function to initialize sigma array (use appropriate logic to populate it)
void initialize_sigma() {
    // This is a simple initialization for demonstration; replace with your actual logic
    for (int mu = 0; mu < MAX_ITER; mu++) {
        for (int deg = 0; deg < MAX_DEG; deg++) {
            for (int alpha_deg = 0; alpha_deg < MAX_ALPHA; alpha_deg++) {
                // Here we're setting some sample values for demonstration.
                // Replace this with your actual computation for sigma coefficients.
                sigma[mu][deg][alpha_deg] = (mu + deg + alpha_deg) % 7;  // Example values
            }
        }
    }
}

// Function to print the sigma coefficients for a given mu
void print_sigma_for_mu(int mu) {
    printf("Sigma coefficients for mu = %d:\n", mu);
    for (int deg = 0; deg < MAX_DEG; deg++) {
        for (int alpha_deg = 0; alpha_deg < MAX_ALPHA; alpha_deg++) {
            uint32_t coefficient = sigma[mu][deg][alpha_deg];
            if (coefficient != 0) {  // We assume 0 means no coefficient at that degree
                printf("sigma[%d][%d][%d] = alpha^%d (value: %d)\n", mu, deg, alpha_deg, coefficient, coefficient);
            }
        }
    }
}

// Function to get the coefficient of x^i in sigma^(mu+1), e.g., for x^3 (i=3)
uint32_t get_coefficient(int mu, int i, int alpha_deg) {
    return sigma[mu][i][alpha_deg];
}

int main() {
    // Initialize the sigma array with sample values
    initialize_sigma();
    
    // Print sigma coefficients for a specific mu (e.g., mu = 2)
    print_sigma_for_mu(2);
    
    // Example of accessing and printing a specific coefficient
    int mu = 2;
    int i = 3;  // Degree of x (for example, x^3)
    int alpha_deg = 2;  // Example alpha degree (e.g., alpha^2)
    uint32_t coefficient = get_coefficient(mu, i, alpha_deg);
    
    printf("Accessed coefficient for sigma[%d][%d][%d] = alpha^%d (value: %d)\n", mu, i, alpha_deg, coefficient, coefficient);
    
    return 0;
}

