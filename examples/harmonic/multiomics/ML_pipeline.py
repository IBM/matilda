from sklearn.preprocessing import StandardScaler, MaxAbsScaler, LabelBinarizer
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    auc,
    roc_curve,
    ConfusionMatrixDisplay,
    confusion_matrix,
    f1_score,
    pairwise_distances,
    accuracy_score,
)
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (
    RandomForestClassifier,
    AdaBoostClassifier,
    GradientBoostingClassifier,
)
from xgboost import XGBClassifier


def ml_pipeline(use_scaler=True):
    # define pipelines for each model
    if use_scaler:
        pipe_nb = Pipeline(steps=[("scaler", StandardScaler()), ("NB", GaussianNB())])
        pipe_knn = Pipeline(
            steps=[("scaler", StandardScaler()), ("KNN", KNeighborsClassifier())]
        )
        pipe_svm = Pipeline(
            steps=[("scaler", StandardScaler()), ("SVM", SVC(random_state=123))]
        )
        pipe_lr = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("LR", LogisticRegression(random_state=123)),
            ]
        )
        pipe_rf = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("RF", RandomForestClassifier(random_state=123)),
            ]
        )
        pipe_ab = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("AB", AdaBoostClassifier(random_state=123)),
            ]
        )
        pipe_grb = Pipeline(
            steps=[
                ("scaler", StandardScaler()),
                ("GBC", GradientBoostingClassifier(random_state=123)),
            ]
        )
        pipe_xgb = Pipeline(
            steps=[("scaler", StandardScaler()), ("XGB", XGBClassifier())]
        )

    else:
        print("no scaler")
        pipe_nb = Pipeline(steps=[("NB", GaussianNB())])
        pipe_knn = Pipeline(steps=[("KNN", KNeighborsClassifier())])
        pipe_svm = Pipeline(steps=[("SVM", SVC(random_state=123))])
        pipe_lr = Pipeline(steps=[("LR", LogisticRegression(random_state=123))])
        pipe_rf = Pipeline(steps=[("RF", RandomForestClassifier(random_state=123))])
        pipe_ab = Pipeline(steps=[("AB", AdaBoostClassifier(random_state=123))])
        pipe_grb = Pipeline(
            steps=[("GBC", GradientBoostingClassifier(random_state=123))]
        )
        pipe_xgb = Pipeline(steps=[("XGB", XGBClassifier())])

    params_nb = [{}]

    params_knn = [
        {
            "KNN__n_neighbors": [2, 3, 4, 5, 6, 7],
            "KNN__weights": ["uniform", "distance"],
            "KNN__metric": ["euclidean", "manhattan", "l1", "l2"],
        }
    ]

    params_svm = [
        {
            "SVM__kernel": ["linear", "rbf", "poly"],
            "SVM__C": [0.1, 0.5, 1, 5, 10],
        }
    ]

    params_lr = [
        {
            "LR__penalty": ["l1", "l2"],
            "LR__C": [0.1, 0.5, 1, 5, 10],
            "LR__solver": ["liblinear"],
        }
    ]

    params_rf = [
        {
            "RF__min_samples_leaf": [1, 2, 3, 4, 5, 6, 7, 8],
            "RF__max_depth": [1, 2, 3, 4, 5, 6, 7],
            "RF__min_samples_split": [2, 3, 4, 5, 6, 7],
            "RF__criterion": ["gini", "entropy"],
        }
    ]

    params_ab = [
        {
            "AB__n_estimators": [10, 20, 40, 50, 70, 100],
            "AB__learning_rate": [0.0001, 0.001, 0.01, 0.1, 1],
        }
    ]

    params_grb = [
        {
            "GBC__n_estimators": [10, 20, 40, 50, 70, 100],
            "GBC__max_depth": range(5, 9, 1),
            "GBC__learning_rate": [0.0001, 0.001, 0.01, 0.1, 1],
            "GBC__criterion": ["friedman_mse", "squared_error"],
            "GBC__min_samples_leaf": [1, 2, 3, 4, 5, 6, 7, 8],
        }
    ]

    params_xgb = [
        {
            # "XGB__n_estimators": [10, 20, 40, 50, 70, 100],
            "XGB__n_estimators": [70, 100],
            "XGB__max_depth": range(5, 9, 1),
            # "XGB__learning_rate": [0.0001, 0.001, 0.01, 0.1, 1],
            "XGB__learning_rate": [0.001, 0.01, 0.1],
            "XGB__min_child_weight": [1, 2, 3, 4, 5],
            "XGB__subsample": [1, 0.5, 0.1],
        }
    ]

    nb_grid_search = GridSearchCV(
        estimator=pipe_nb,
        param_grid=params_nb,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    knn_grid_search = GridSearchCV(
        estimator=pipe_knn,
        param_grid=params_knn,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    svm_grid_search = GridSearchCV(
        estimator=pipe_svm,
        param_grid=params_svm,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    lr_grid_search = GridSearchCV(
        estimator=pipe_lr,
        param_grid=params_lr,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    rf_grid_search = GridSearchCV(
        estimator=pipe_rf,
        param_grid=params_rf,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    ab_grid_search = GridSearchCV(
        estimator=pipe_ab,
        param_grid=params_ab,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    grb_grid_search = GridSearchCV(
        estimator=pipe_grb,
        param_grid=params_grb,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    xgb_grid_search = GridSearchCV(
        estimator=pipe_xgb,
        param_grid=params_xgb,
        scoring="f1_micro",
        cv=5,
        verbose=1,
        n_jobs=-1,
    )

    all_grids = [
        # nb_grid_search,
        # knn_grid_search,
        # svm_grid_search,
        # lr_grid_search,
        # rf_grid_search,
        # ab_grid_search,
        # grb_grid_search,
        xgb_grid_search,
    ]

    return all_grids
