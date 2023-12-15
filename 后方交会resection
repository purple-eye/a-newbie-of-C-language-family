from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QPushButton,
    QLabel,
    QLineEdit,
    QVBoxLayout,
    QHBoxLayout,
)
import numpy as np
import sys


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Resection")
        self.setGeometry(0, 0, 800, 600)
        topLayout = QVBoxLayout(self)
        self.setLayout(topLayout)
        inLayout = QVBoxLayout()
        topLayout.addLayout(inLayout)

        h1 = QHBoxLayout()
        inLayout.addLayout(h1)

        self.inputx0 = QLineEdit()
        h1.addWidget(QLabel("x0:"))
        h1.addWidget(self.inputx0)

        self.inputy0 = QLineEdit()
        h2 = QHBoxLayout()
        inLayout.addLayout(h2)
        h2.addWidget(QLabel("y0:"))
        h2.addWidget(self.inputy0)

        self.inputf = QLineEdit()
        h3 = QHBoxLayout()
        inLayout.addLayout(h3)
        h3.addWidget(QLabel("f:"))
        h3.addWidget(self.inputf)

        self.m = QLineEdit()
        h4 = QHBoxLayout()
        inLayout.addLayout(h4)
        h4.addWidget(QLabel("m:"))
        h4.addWidget(self.m)

        coordinateLayout = QHBoxLayout()
        topLayout.addLayout(coordinateLayout)

        xLayout = QVBoxLayout()
        coordinateLayout.addLayout(xLayout)
        xLayout.addWidget(QLabel("x(mm)"))
        self.inputx1 = QLineEdit()
        xLayout.addWidget(self.inputx1)
        self.inputx2 = QLineEdit()
        xLayout.addWidget(self.inputx2)
        self.inputx3 = QLineEdit()
        xLayout.addWidget(self.inputx3)
        self.inputx4 = QLineEdit()
        xLayout.addWidget(self.inputx4)

        yLayout = QVBoxLayout()
        coordinateLayout.addLayout(yLayout)
        yLayout.addWidget(QLabel("y(mm)"))
        self.inputy1 = QLineEdit()
        yLayout.addWidget(self.inputy1)
        self.inputy2 = QLineEdit()
        yLayout.addWidget(self.inputy2)
        self.inputy3 = QLineEdit()
        yLayout.addWidget(self.inputy3)
        self.inputy4 = QLineEdit()
        yLayout.addWidget(self.inputy4)
        XLayout = QVBoxLayout()
        coordinateLayout.addLayout(XLayout)
        XLayout.addWidget(QLabel("X(m)"))
        self.inputX1 = QLineEdit()
        XLayout.addWidget(self.inputX1)
        self.inputX2 = QLineEdit()
        XLayout.addWidget(self.inputX2)
        self.inputX3 = QLineEdit()
        XLayout.addWidget(self.inputX3)
        self.inputX4 = QLineEdit()
        XLayout.addWidget(self.inputX4)

        YLayout = QVBoxLayout()
        coordinateLayout.addLayout(YLayout)
        YLayout.addWidget(QLabel("Y(m)"))
        self.inputY1 = QLineEdit()
        YLayout.addWidget(self.inputY1)
        self.inputY2 = QLineEdit()
        YLayout.addWidget(self.inputY2)
        self.inputY3 = QLineEdit()
        YLayout.addWidget(self.inputY3)
        self.inputY4 = QLineEdit()
        YLayout.addWidget(self.inputY4)

        ZLayout = QVBoxLayout()
        coordinateLayout.addLayout(ZLayout)
        ZLayout.addWidget(QLabel("Z(m)"))
        self.inputZ1 = QLineEdit()
        ZLayout.addWidget(self.inputZ1)
        self.inputZ2 = QLineEdit()
        ZLayout.addWidget(self.inputZ2)
        self.inputZ3 = QLineEdit()
        ZLayout.addWidget(self.inputZ3)
        self.inputZ4 = QLineEdit()
        ZLayout.addWidget(self.inputZ4)
        self.resultlabel = QLabel("", self)
        topLayout.addWidget(self.resultlabel)

        button = QPushButton("开始计算")
        button.clicked.connect(self.c)
        topLayout.addWidget(button)

    def c(self):
        # 输入的数据
        x0 = float(self.inputx0.text())
        y0 = float(self.inputy0.text())
        f = float(self.inputf.text())
        x = [
            float(self.inputx1.text()),
            float(self.inputx2.text()),
            float(self.inputx3.text()),
            float(self.inputx4.text()),
        ]
        y = [
            float(self.inputy1.text()),
            float(self.inputy2.text()),
            float(self.inputy3.text()),
            float(self.inputy4.text()),
        ]
        X = [
            float(self.inputX1.text()),
            float(self.inputX2.text()),
            float(self.inputX3.text()),
            float(self.inputX4.text()),
        ]
        Y = [
            float(self.inputY1.text()),
            float(self.inputY2.text()),
            float(self.inputY3.text()),
            float(self.inputY4.text()),
        ]
        Z = [
            float(self.inputZ1.text()),
            float(self.inputZ2.text()),
            float(self.inputZ3.text()),
            float(self.inputZ4.text()),
        ]
        m = float(self.m.text())
        x = np.array(x)
        y = np.array(y)
        X = np.array(X)
        Y = np.array(Y)

        # 确定未知数的初值
        self.Xs0 = X.mean()
        self.Ys0 = Y.mean()
        self.Zs0 = m * f
        self.f0 = 0
        self.w0 = 0
        self.k0 = 0
        x_apxm = [0, 0, 0, 0]
        y_apxm = [0, 0, 0, 0]

        R = np.mat(np.zeros((3, 3)))
        L = np.mat(np.zeros((8, 1)))
        A = np.mat(np.zeros((8, 6)))
        result = ""
        f_cor = 1
        w_cor = 1
        k_cor = 1
        flag = 0
        # 迭代计算
        while flag < 100:
            flag += 1
            if (
                (abs(f_cor) > 0.000001)
                | (abs(w_cor) > 0.000001)
                | (abs(k_cor) > 0.000001)
            ):
                R = self.r_mat()
                x_apxm, y_apxm = self.xy_approximate(X, Y, Z, x, y, R, x0, y0, f)

                for i in range(4):
                    L[2 * i] = x_apxm[i]
                    L[2 * i + 1] = y_apxm[i]

                for i in range(4):
                    A[2 * i : 2 * i + 2, :] = self.a_parameter(
                        X[i],
                        Y[i],
                        Z[i],
                        self.Xs0,
                        self.Ys0,
                        self.Zs0,
                        x[i],
                        y[i],
                        self.w0,
                        self.k0,
                        R,
                        f,
                        x0,
                        y0,
                    )

                X_mat = np.mat(np.zeros((6, 1)))
                X_mat = (A.T * A).I * A.T * L

                f_cor = X_mat[3, 0]
                w_cor = X_mat[4, 0]
                k_cor = X_mat[5, 0]

                self.Xs0 = self.Xs0 + X_mat[0, 0]
                self.Ys0 = self.Ys0 + X_mat[1, 0]
                self.Zs0 = self.Zs0 + X_mat[2, 0]
                self.f0 = self.f0 + X_mat[3, 0]
                self.w0 = self.w0 + X_mat[4, 0]
                self.k0 = self.k0 + X_mat[5, 0]
                V = np.mat(np.zeros((8, 1)))
                V = A * X_mat - L
                errorValue = np.sqrt((V.T * V) / (2 * 4 - 6))
                result = "最终结果为：Xs=%f, Ys=%f, Zs=%f, f=%f, w=%f, k=%f,\n单位权中误差为%f " % (
                    self.Xs0,
                    self.Ys0,
                    self.Zs0,
                    self.f0,
                    self.w0,
                    self.k0,
                    errorValue,
                )
            else:
                self.resultlabel.setText(result)
                break

    # 计算旋转矩阵
    def r_mat(self):
        Rf = np.mat(
            [
                [np.cos(self.f0), 0, -np.sin(self.f0)],
                [0.0, 1, 0],
                [np.sin(self.f0), 0, np.cos(self.f0)],
            ]
        )
        Rw = np.mat(
            [
                [1.0, 0, 0],
                [0, np.cos(self.w0), -np.sin(self.w0)],
                [0, np.sin(self.w0), np.cos(self.w0)],
            ]
        )
        Rk = np.mat(
            [
                [np.cos(self.k0), -np.sin(self.k0), 0],
                [np.sin(self.k0), np.cos(self.k0), 0],
                [0, 0.0, 1],
            ]
        )
        R = Rf * Rw * Rk
        return R

    # 计算x，y近似值
    def xy_approximate(self, X, Y, Z, x, y, R, x0, y0, f):
        self.x_apxm = [0, 0, 0, 0]
        self.y_apxm = [0, 0, 0, 0]

        for i in range(4):
            self.x_apxm[i] = x[i] - (
                x0
                - f
                * (
                    (
                        R[0, 0] * (X[i] - self.Xs0)
                        + R[1, 0] * (Y[i] - self.Ys0)
                        + R[2, 0] * (Z[i] - self.Zs0)
                    )
                    / (
                        R[0, 2] * (X[i] - self.Xs0)
                        + R[1, 2] * (Y[i] - self.Ys0)
                        + R[2, 2] * (Z[i] - self.Zs0)
                    )
                )
            )
            self.y_apxm[i] = y[i] - (
                y0
                - f
                * (
                    (
                        R[0, 1] * (X[i] - self.Xs0)
                        + R[1, 1] * (Y[i] - self.Ys0)
                        + R[2, 1] * (Z[i] - self.Zs0)
                    )
                    / (
                        R[0, 2] * (X[i] - self.Xs0)
                        + R[1, 2] * (Y[i] - self.Ys0)
                        + R[2, 2] * (Z[i] - self.Zs0)
                    )
                )
            )

        return self.x_apxm, self.y_apxm

    # 计算系数矩阵
    def a_parameter(self, X, Y, Z, Xs, Ys, Zs, x, y, w, k, R, f, x0, y0):
        self.parameter = np.zeros((2, 6))
        mean = np.zeros((3, 1))
        minus = np.zeros((3, 1))

        minus = np.array([[X - Xs], [Y - Ys], [Z - Zs]])
        mean = R.T * np.mat(minus)

        self.parameter[0][0] = (R[0, 0] * f + R[0, 2] * (x - x0)) / mean[2]
        self.parameter[0][1] = (R[1, 0] * f + R[1, 2] * (x - x0)) / mean[2]
        self.parameter[0][2] = (R[2, 0] * f + R[2, 2] * (x - x0)) / mean[2]
        self.parameter[1][0] = (R[0, 1] * f + R[0, 2] * (y - y0)) / mean[2]
        self.parameter[1][1] = (R[1, 1] * f + R[1, 2] * (y - y0)) / mean[2]
        self.parameter[1][2] = (R[2, 1] * f + R[2, 2] * (y - y0)) / mean[2]

        self.parameter[0][3] = (y - y0) * np.sin(w) - (
            ((x - x0) / f) * ((x - x0) * np.cos(k) - (y - y0) * np.sin(k))
            + f * np.cos(k)
        ) * np.cos(w)
        self.parameter[0][4] = -f * np.sin(k) - ((x - x0) / f) * (
            (x - x0) * np.sin(k) + (y - y0) * np.cos(k)
        )
        self.parameter[0][5] = y - y0
        self.parameter[1][3] = -(x - x0) * np.sin(w) - (
            ((y - y0) / f) * ((x - x0) * np.cos(k) - (y - y0) * np.sin(k))
            - f * np.cos(k)
        ) * np.cos(w)
        self.parameter[1][4] = -f * np.cos(k) - ((y - y0) / f) * (
            (x - x0) * np.sin(k) + (y - y0) * np.cos(k)
        )
        self.parameter[1][5] = -(x - x0)

        return self.parameter


app = QApplication(sys.argv)
mw = MainWindow()
mw.show()
sys.exit(app.exec_())
