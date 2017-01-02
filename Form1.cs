using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

namespace SOM
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            textBox1.Clear();
            Bitmap currentPicture = new Bitmap(System.Drawing.Image.FromFile("0.png", true));
            int inputNeurons = currentPicture.Size.Height * currentPicture.Size.Width;
            int outputNeurons = 6;
            double[][] InputVector = new double[outputNeurons][];
            for (int i = 0; i < outputNeurons; i++)
            {
                currentPicture = new Bitmap(System.Drawing.Image.FromFile(AppDomain.CurrentDomain.BaseDirectory + i.ToString() + ".png", true));
                InputVector[i] = ImageToVector(currentPicture);
            }
            NeuroNet net = new NeuroNet(inputNeurons, outputNeurons);
            Study(ref net, InputVector);
            for (int i = 0; i < 6; i++)
            {
                textBox1.Text += "При подаче на вход " + i + "-го образца, Ответ: " + Test(net, InputVector[i]).ToString() + Environment.NewLine;
            }
        }

        private double[] ImageToVector(Bitmap img)
        {
            int Size = img.Size.Height * img.Size.Width;
            double[] vector = new double[Size];
            int i = 0;
            for (int x = 0; x < img.Size.Width; x++)
            {
                for (int y = 0; y < img.Size.Height; y++)
                {
                    Color pixel = img.GetPixel(x, y);
                    Byte lum = (Byte)((pixel.R * 77 + pixel.G * 151 + pixel.B * 28) >> 8);
                    vector[i++] = 1.0f - lum / 255.0f;
                }
            }
            return vector;
        }

        private int Test(NeuroNet net, double[] InputVector)
        {
            double MinDistance = EuclideanDistance(net.neurons[0], InputVector);
            int BMUIndex = 0;
            for (int i = 1; i < net.neurons.Count; i++)
            {
                double tmp_ED = EuclideanDistance(net.neurons[i], InputVector);
                if (tmp_ED < MinDistance)
                {
                    BMUIndex = i;
                    MinDistance = tmp_ED;
                }
            }
            return BMUIndex;
        }

        private void Study(ref NeuroNet net, double[][] InputVector)
        {
            int c;
            for (int k = 0; k < 6; k++) // цикл, в котором предъявляем сети входные вектора - InputVector
            {
                double MinDistance = EuclideanDistance(net.neurons[0], InputVector[k]);
                int BMUIndex = 0;
                for (int i = 1; i < net.neurons.Count; i++)
                {
                    double tmp_ED = EuclideanDistance(net.neurons[i], InputVector[k]); // Находим Евклидово расстояние между i-ым нейроном и k-ым входным  вектором
                    if (tmp_ED < MinDistance) // Eсли Евклидово расстояние минимально, то это нейрон-победитель
                    {
                        BMUIndex = i; // Индекс нейрона-победителя
                        MinDistance = tmp_ED;
                    }
                }

                for (int i = 0; i < net.neurons.Count; i++)
                {
                    for (int g = 0; g < InputVector[k].Length; g++)
                    {
                        double hfunc = hc(k, net.neurons[BMUIndex].weights[g], net.neurons[i].weights[g]);
                        double normfunc = normLearningRate(k);
                        //net.neurons[i].weights[g] = net.neurons[i].weights[g] + hfunc * normfunc * (InputVector[k][g] - net.neurons[i].weights[g]);
                        net.neurons[i].weights[g] = net.neurons[i].weights[g] + hfunc * normfunc * (InputVector[i][g] - net.neurons[i].weights[g]); // Вроде так работает
                        if (i > 0 && g > 282)
                            c = 0;
                    }
                }
                double Error = EuclideanDistance(net.neurons[BMUIndex], InputVector[k]);
                for (int y = 0; y < 6; y++)
                {
                    textBox1.Text += "Евклидово расстояние " + y + "-го нейрона между " + y + "-м входным вектором на " + k + "-ой итерации: " + EuclideanDistance(net.neurons[y], InputVector[y]) + Environment.NewLine;
                }
                textBox1.Text += Environment.NewLine;
            }
        }

        private double hc(int k, double winnerCoordinate, double Coordinate)
        {
            double dist = Distance(winnerCoordinate, Coordinate);
            //double s = sigma(k);
            return Math.Exp(-dist / 2 * Sqr(sigma(k)));
        }

        private double sigma(int k)
        {
            //return -0.01 * k + 2;
            return 1 * Math.Exp(-k / 5);

            //double nf = 1000 / Math.Log(2025);
            //return Math.Exp(-k / nf) * 2025;
        }

        private double normLearningRate(int k)
        {
            return 0.1 * Math.Exp(-k / 1000);
        }

        private double EuclideanDistance(Neuron neuron, double[] InputVector)
        {
            double Sum = 0;
            for (int i = 0; i < InputVector.Length; i++)
            {
                Sum += Sqr(InputVector[i] - neuron.weights[i]);
            }
            return Math.Sqrt(Sum);
        }

        private double Distance(double winnerCoordinate, double Coordinate)
        {
            return Math.Sqrt(Sqr(winnerCoordinate - Coordinate));
        }

        private double Sqr(double value)
        {
            return value * value;
        }
    }

    public class NeuroNet
    {
        public int inputs = 0;
        public List<Neuron> neurons;

        public NeuroNet(int inputs_, int neurons_)
        {
            neurons = new List<Neuron>();
            inputs = inputs_;
            for (int i = 0; i < neurons_; i++)
            {
                Neuron neuron = new Neuron(i, inputs_);
                neurons.Add(neuron);
            }
        }
    }


    public class Neuron
    {
        public int number = 0;
        public List<double> weights;

        public Neuron(int number_, int inputs_)
        {
            weights = new List<double>();
            number = number_;
            Random rand = new Random();
            for (int i = 0; i < inputs_; i++)
            {
                weights.Add(rand.NextDouble());
            }

        }
    }
}