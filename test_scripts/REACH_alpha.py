#!/usr/bin/env python3

"""
Enhancer Prediction Pipeline
This script identifies enhancer elements in DNA using supervised learning models
and multiple data types including Hi-C, ChIP-seq, and RNA-seq data.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, cross_val_score, GridSearchCV
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from sklearn.metrics import classification_report, roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import StandardScaler
import pyBigWig
import cooler
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional
import requests
import logging
import sys
import json
from pathlib import Path
import urllib.parse
import time
from concurrent.futures import ThreadPoolExecutor
import warnings
warnings.filterwarnings('ignore')

class EnhancerVisualizer:
    """Class to handle all visualization tasks"""
    
    def __init__(self, output_dir: str = 'enhancer_analysis'):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
    
    def plot_genomic_distribution(self, enhancers: pd.DataFrame):
        """Plot genomic distribution of enhancers"""
        plt.figure(figsize=(15, 5))
        
        plt.subplot(131)
        enhancers['chromosome'].value_counts().plot(kind='bar')
        plt.title('Enhancer Distribution by Chromosome')
        plt.xlabel('Chromosome')
        plt.ylabel('Count')
        
        plt.subplot(132)
        enhancer_lengths = enhancers['end'] - enhancers['start']
        sns.histplot(enhancer_lengths)
        plt.title('Enhancer Length Distribution')
        plt.xlabel('Length (bp)')
        
        plt.subplot(133)
        enhancers['source'].value_counts().plot(kind='pie', autopct='%1.1f%%')
        plt.title('Enhancer Source Distribution')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'genomic_distribution.png')
        plt.close()
    
    def plot_feature_distributions(self, X: np.ndarray, feature_names: List[str]):
        """Plot distribution of features"""
        n_features = X.shape[1]
        n_cols = 3
        n_rows = (n_features + n_cols - 1) // n_cols
        
        plt.figure(figsize=(15, 5 * n_rows))
        
        for i, feature in enumerate(feature_names):
            plt.subplot(n_rows, n_cols, i + 1)
            sns.histplot(X[:, i], kde=True)
            plt.title(f'{feature} Distribution')
            
        plt.tight_layout()
        plt.savefig(self.output_dir / 'feature_distributions.png')
        plt.close()
    
    def plot_interaction_heatmap(self, interaction_matrix: np.ndarray, 
                               region: tuple, target: tuple):
        """Plot Hi-C interaction heatmap"""
        plt.figure(figsize=(8, 8))
        sns.heatmap(interaction_matrix, cmap='YlOrRd')
        plt.title(f'Hi-C Interactions\n{region[0]}:{region[1]}-{region[2]} vs {target[0]}:{target[1]}-{target[2]}')
        plt.xlabel('Target Region Bins')
        plt.ylabel('Enhancer Region Bins')
        plt.savefig(self.output_dir / 'interaction_heatmap.png')
        plt.close()
    
    def plot_model_comparison(self, results: pd.DataFrame):
        """Plot model comparison results"""
        plt.figure(figsize=(12, 6))
        
        plt.subplot(121)
        sns.barplot(data=results, x='Model', y='Mean CV Score')
        plt.xticks(rotation=45)
        plt.title('Model Performance Comparison')
        
        plt.subplot(122)
        sns.barplot(data=results, x='Model', y='Training Time')
        plt.xticks(rotation=45)
        plt.title('Model Training Time Comparison')
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'model_comparison.png')
        plt.close()

class EnhancerPredictor:
    """Main class for enhancer prediction"""
    
    def __init__(self, hic_file: str, chip_files: Dict[str, str], 
                 rna_file: str, genome_file: str):
        self._setup_logging()
        try:
            self.hic = cooler.Cooler(hic_file)
            self.chip_data = {mark: pyBigWig.open(path) for mark, path in chip_files.items()}
            self.expression = pd.read_csv(rna_file, index_col=0)
            self.genome = pd.read_csv(genome_file, sep='\t')
            self.model = None
            self.feature_names = None
            self.scaler = StandardScaler()
            self.visualizer = EnhancerVisualizer()
            self.data_fetcher = EnhancerDataFetcher()
            self.model_evaluator = ModelEvaluator()
        except Exception as e:
            self.logger.error(f"Failed to initialize: {str(e)}")
            raise

    def _setup_logging(self):
        """Configure logging"""
        self.logger = logging.getLogger(__name__)
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        self.logger.addHandler(handler)
        self.logger.setLevel(logging.INFO)



    def create_training_data(self, enhancers: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Create training dataset"""
        features = []
        labels = []
        
        for _, row in enhancers.iterrows():
            region = (row['chromosome'], row['start'], row['end'])
            promoter = self.get_promoter_region(row['target_gene'])
            
            hic_feats = self.extract_hic_features(region, promoter)
            chip_feats = self.extract_chip_features(region)
            expr_change = self.get_expression_change(row['target_gene'], 'sample1', 'sample2')
            
            combined_feats = np.concatenate([hic_feats, chip_feats, [expr_change]])
            features.append(combined_feats)
            labels.append(1 if row['confidence'] == 'Validated' else 0)
        
        return np.array(features), np.array(labels)

def run_complete_analysis(self, output_dir: str = 'enhancer_analysis'):
    """
    Run a comprehensive enhancer prediction and analysis pipeline
    
    This method performs the following steps:
    1. Fetch enhancers from multiple databases
    2. Visualize genomic distribution
    3. Create training dataset
    4. Evaluate and compare machine learning models
    5. Select and save the best model
    6. Generate detailed visualizations
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    # Log start of analysis
    self.logger.info("Starting comprehensive enhancer analysis pipeline")

    try:
        # 1. Fetch enhancers from multiple databases
        enhancers_vista = self.data_fetcher.fetch_vista_enhancers()
        enhancers_encode = self.data_fetcher.fetch_encode_enhancers()
        enhancers_genehancer = self.data_fetcher.fetch_genehancer_enhancers()

        # Combine enhancers from different sources
        all_enhancers = pd.concat([
            enhancers_vista, 
            enhancers_encode, 
            enhancers_genehancer
        ], ignore_index=True)

        # 2. Visualize genomic distribution
        self.visualizer.plot_genomic_distribution(all_enhancers)

        # 3. Create training dataset
        X, y = self.create_training_data(all_enhancers)

        # Extract feature names for visualization and logging
        self.feature_names = [
            # Hi-C features
            'hic_mean', 'hic_max', 'hic_sum', 'hic_std',
            # ChIP-seq features (repeated for each mark)
            'h3k4me1_mean', 'h3k4me1_max', 'h3k4me1_std',
            'h3k27ac_mean', 'h3k27ac_max', 'h3k27ac_std',
            # Add more marks as needed
            'expression_change'
        ]

        # 4. Visualize feature distributions
        self.visualizer.plot_feature_distributions(X, self.feature_names)

        # 5. Scale features
        X_scaled = self.scaler.fit_transform(X)

        # 6. Split data into training and testing sets
        X_train, X_test, y_train, y_test = train_test_split(
            X_scaled, y, test_size=0.2, random_state=42
        )

        # 7. Evaluate and compare models
        model_results = self.model_evaluator.evaluate_all_models(X_train, y_train)

        # 8. Visualize model comparison
        self.visualizer.plot_model_comparison(model_results)

        # 9. Select the best model
        best_model_row = model_results.loc[model_results['Mean CV Score'].idxmax()]
        self.model = best_model_row['Fitted Model']
        
        # 10. Evaluate on test set
        y_pred = self.model.predict(X_test)
        y_pred_proba = self.model.predict_proba(X_test)[:, 1]

        # 11. Generate classification report
        report = classification_report(y_test, y_pred)
        with open(output_path / 'classification_report.txt', 'w') as f:
            f.write(report)

        # 12. Plot ROC curve
        fpr, tpr, _ = roc_curve(y_test, y_pred_proba)
        roc_auc = auc(fpr, tpr)

        plt.figure()
        plt.plot(fpr, tpr, color='darkorange', lw=2, 
                 label=f'ROC curve (AUC = {roc_auc:.2f})')
        plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver Operating Characteristic (ROC)')
        plt.legend(loc="lower right")
        plt.savefig(output_path / 'roc_curve.png')
        plt.close()

        # 13. Save model and scaler
        import joblib
        joblib.dump(self.model, output_path / 'best_enhancer_model.pkl')
        joblib.dump(self.scaler, output_path / 'feature_scaler.pkl')

        self.logger.info("Enhancer analysis pipeline completed successfully")
        return model_results, report

    except Exception as e:
        self.logger.error(f"Error in analysis pipeline: {str(e)}")
        raise
