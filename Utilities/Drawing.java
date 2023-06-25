package Utilities;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import javax.imageio.ImageIO;
import javax.swing.*;

public class Drawing {

    /**
     * Enum representing the type of a 3D rectangle (inward or outward).
     */
    public enum Type {
        in, out
    }

    private JFrame frame = new JFrame();

    /**
     * Default constructor for the Drawing class. Creates a JFrame with default
     * settings.
     */
    public Drawing() {
        configureJFrame(400, 400);
    }

    /**
     * Constructor for the Drawing class with custom dimensions. Creates a
     * JFrame with the specified width and height.
     *
     * @param width The width of the JFrame.
     * @param height The height of the JFrame.
     */
    public Drawing(int width, int height) {
        configureJFrame(width, height);
    }

    /**
     * Configures the JFrame with the specified width and height.
     *
     * @param width The width of the JFrame.
     * @param height The height of the JFrame.
     */
    private void configureJFrame(int width, int height) {
        frame.setVisible(true);
        frame.setResizable(false);
        frame.setSize(width, height);
        frame.setLocationRelativeTo(null);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }

    /**
     * Adds a panel to the JFrame.
     *
     * @param panel The panel to be added.
     */
    private void addPanelToFrame(JPanel panel) {
        frame.add(panel);
        panel.setBounds(0, 0, frame.getWidth(), frame.getHeight());
    }

    /**
     * Draws a line on the JFrame.
     *
     * @param border The width of the line.
     * @param startX The x-coordinate of the starting point.
     * @param startY The y-coordinate of the starting point.
     * @param endX The x-coordinate of the ending point.
     * @param endY The y-coordinate of the ending point.
     * @param color The color of the line.
     */
    public void drawLine(int border, int startX, int startY, int endX, int endY, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawLine(startX, startY, endX, endY);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a 3D rectangle on the JFrame.
     *
     * @param startX The x-coordinate of the upper-left corner of the rectangle.
     * @param startY The y-coordinate of the upper-left corner of the rectangle.
     * @param endX The width of the rectangle.
     * @param endY The height of the rectangle.
     * @param color The color of the rectangle.
     * @param type The type of the rectangle (inward or outward).
     */
    public void draw3DRect(int startX, int startY, int endX, int endY, Color color, Type type) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setColor(color);
                boolean isOutward = (type != Type.in);
                g2.draw3DRect(startX, startY, endX, endY, isOutward);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws an arc on the JFrame.
     *
     * @param border The width of the arc.
     * @param startX The x-coordinate of the upper-left corner of the arc.
     * @param startY The y-coordinate of the upper-left corner of the arc.
     * @param endX The width of the arc.
     * @param endY The height of the arc.
     * @param startTheta The starting angle of the arc in degrees.
     * @param theta The extent of the arc in degrees.
     * @param color The color of the arc.
     */
    public void drawArc(int border, int startX, int startY, int endX, int endY, int startTheta, int theta, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawArc(startX, startY, endX, endY, startTheta, theta);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws an oval on the JFrame.
     *
     * @param border The width of the oval.
     * @param startX The x-coordinate of the upper-left corner of the oval.
     * @param startY The y-coordinate of the upper-left corner of the oval.
     * @param width The width of the oval.
     * @param height The height of the oval.
     * @param color The color of the oval.
     */
    public void drawOval(int border, int startX, int startY, int width, int height, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawOval(startX, startY, width, height);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a polygon on the JFrame.
     *
     * @param border The width of the polygon lines.
     * @param xPoints The array of x-coordinates of the polygon vertices.
     * @param yPoints The array of y-coordinates of the polygon vertices.
     * @param nPoints The number of polygon vertices.
     * @param color The color of the polygon.
     */
    public void drawPolygon(int border, int[] xPoints, int[] yPoints, int nPoints, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawPolygon(xPoints, yPoints, nPoints);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a polyline on the JFrame.
     *
     * @param border The width of the polyline lines.
     * @param xPoints The array of x-coordinates of the polyline vertices.
     * @param yPoints The array of y-coordinates of the polyline vertices.
     * @param nPoints The number of polyline vertices.
     * @param color The color of the polyline.
     */
    public void drawPolyline(int border, int[] xPoints, int[] yPoints, int nPoints, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawPolyline(xPoints, yPoints, nPoints);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a rectangle on the JFrame.
     *
     * @param border The width of the rectangle lines.
     * @param startX The x-coordinate of the upper-left corner of the rectangle.
     * @param startY The y-coordinate of the upper-left corner of the rectangle.
     * @param endX The width of the rectangle.
     * @param endY The height of the rectangle.
     * @param color The color of the rectangle.
     */
    public void drawRect(int border, int startX, int startY, int endX, int endY, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawRect(startX, startY, endX, endY);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a rounded rectangle on the JFrame.
     *
     * @param border The width of the rounded rectangle lines.
     * @param startX The x-coordinate of the upper-left corner of the rounded
     * rectangle.
     * @param startY The y-coordinate of the upper-left corner of the rounded
     * rectangle.
     * @param width The width of the rounded rectangle.
     * @param height The height of the rounded rectangle.
     * @param radius The corner radius of the rounded rectangle.
     * @param color The color of the rounded rectangle.
     */
    public void drawRoundRect(int border, int startX, int startY, int width, int height, int radius, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setStroke(new BasicStroke(border));
                g2.setColor(color);
                g2.drawRoundRect(startX, startY, width, height, radius, radius);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a string on the JFrame.
     *
     * @param fontSize The size of the font.
     * @param text The string to be drawn.
     * @param x The x-coordinate of the starting point of the string.
     * @param y The y-coordinate of the starting point of the string.
     * @param color The color of the string.
     */
    public void drawString(int fontSize, String text, int x, int y, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                g2.setFont(new Font(Font.SANS_SERIF, Font.PLAIN, fontSize));
                g2.setColor(color);
                g2.drawString(text, x, y);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Draws a clear rectangle on the JFrame.
     *
     * @param startX The x-coordinate of the upper-left corner of the rectangle.
     * @param startY The y-coordinate of the upper-left corner of the rectangle.
     * @param width The width of the rectangle.
     * @param height The height of the rectangle.
     * @param color The background color of the rectangle.
     */
    public void drawClearRect(int startX, int startY, int width, int height, Color color) {
        JPanel Panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                g.setColor(color);
                g.clearRect(startX, startY, width, height);
            }
        };
        addPanelToFrame(Panel);
    }

    /**
     * Saves the drawing as an image file.
     *
     * @param filePath The file path to save the image.
     * @param format The image format (e.g., "png", "jpg").
     */
    public void saveDrawing(String filePath, String format) {
        try {
            Container content = frame.getContentPane();
            BufferedImage image = new BufferedImage(content.getWidth(), content.getHeight(), BufferedImage.TYPE_INT_RGB);
            Graphics2D g2d = image.createGraphics();
            content.printAll(g2d);
            g2d.dispose();
            ImageIO.write(image, format, new File(filePath));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Draws an image on a JPanel at the specified coordinates.
     *
     * @param startX The x-coordinate to start drawing the image.
     * @param startY The y-coordinate to start drawing the image.
     * @param URL The file path or URL of the image.
     */
    public void drawImage(int startX, int startY, String URL) {
        JPanel panel = new JPanel() {
            @Override
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;
                File imageFile = new File(URL);
                try {
                    Image image = ImageIO.read(imageFile);
                    g2.drawImage(image, startX, startY, this);
                } catch (Exception ex) {
                    JOptionPane.showMessageDialog(null, ex);
                }
            }
        };
        frame.add(panel);
        panel.setBounds(0, 0, frame.getWidth(), frame.getHeight());
    }

    public void drawCoordinateLines(int width) {

        //currently not working
        int mid = width / 2;

        // Draw x-axis
        drawLine(2, 0, mid, width, mid, Color.BLACK);

        // Draw y-axis
        drawLine(2, mid, 0, mid, width, Color.BLACK);

        int labelOffset = 15;
        int tickSize = 1;

        // Draw positive y-axis (+y)
        for (int y = mid, i = 1; y >= 0; y -= 25, i++) {
            drawLine(0, y, width, y, tickSize, Color.CYAN);
            drawString(15, Integer.toString(i), mid - labelOffset, y + labelOffset, Color.BLACK);
        }

        // Draw negative y-axis (-y)
        for (int y = mid, i = -1; y < width; y += 25, i--) {
            drawLine(0, y, width, y, tickSize, Color.CYAN);
            drawString(15, Integer.toString(i), mid - labelOffset, y + labelOffset, Color.BLACK);
        }

        // Draw negative x-axis (x-)
        for (int x = mid, i = -1; x >= 0; x -= 25, i--) {
            drawLine(x, 0, x, width, tickSize, Color.CYAN);
            drawString(15, Integer.toString(i), x - labelOffset, mid + labelOffset, Color.BLACK);
        }

        // Draw positive x-axis (x+)
        for (int x = mid, i = 1; x < width; x += 25, i++) {
            drawLine(x, 0, x, width, tickSize, Color.CYAN);
            drawString(15, Integer.toString(i), x - labelOffset, mid + labelOffset, Color.BLACK);
        }
    }

    /**
     * Draws a linear function on the graph.
     *
     * @param width The width of the graph.
     * @param a The coefficient of x in the linear function.
     * @param c The constant term in the linear function.
     */
    public void drawLinearFunction(int width, int a, int c) {
        frame.setSize(width, width);
        int nPoints = frame.getWidth();
        int mid = nPoints / 2;
        int[] xPoints = new int[nPoints];
        int[] yPoints = new int[nPoints];

        for (int x = 0, i = -mid; x < nPoints; x++, i++) {
            xPoints[x] = x;
            double fx = a * i + c * 25;
            yPoints[x] = (int) -fx + mid;
        }

        drawPolyline(0, xPoints, yPoints, nPoints, Color.BLUE);
        drawCoordinateLines(width);
    }

    /**
     * Draws a quadratic function on the graph.
     *
     * @param width The width of the graph.
     * @param a The coefficient of x^2 in the quadratic function.
     * @param b The coefficient of x in the quadratic function.
     * @param c The constant term in the quadratic function.
     * @param points The number of points to use for plotting the function.
     */
    public void drawQuadraticFunction(int width, int a, int b, int c, int points) {
        frame.setSize(width, width);
        frame.setPreferredSize(new Dimension(width, width));
        int[] xPoints = new int[points];
        int[] yPoints = new int[points];
        int mid = width / 2;
        double step = (double) width / (points - 1);

        for (int i = 0; i < points; i++) {
            double x = i * step;
            xPoints[i] = (int) x;
            double fx = a * Math.pow((x - mid) / 25, 2) + b * ((x - mid) / 25) + c;
            yPoints[i] = (int) (-fx * 25) + mid;
        }
        drawCoordinateLines(width);
        drawPolyline(0, xPoints, yPoints, points, Color.BLUE);

    }

}
