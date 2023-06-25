package CustomGUI;

import java.awt.Color;
import java.awt.Cursor;
import java.awt.Graphics;

import javax.swing.JComboBox;

/**
 * A custom combo box component with a rounded rectangle shape and a blue
 * border.
 */
public class CustomComboBox extends JComboBox {

    private int R;
    private int B;
    private int G;
    private Color border;

    /**
     * Constructs a new CustomComboBox instance with default values.
     */
    public CustomComboBox() {
        this(new Cursor(Cursor.HAND_CURSOR), 70, 130, 180, Color.BLUE);
    }

    /**
     * Constructs a new CustomComboBox instance with the specified parameters.
     *
     * @param cursor the cursor to be set for the combo box
     * @param R the red component value for the background color
     * @param B the blue component value for the background color
     * @param G the green component value for the background color
     * @param border the color for the border
     */
    public CustomComboBox(Cursor cursor, int R, int B, int G, Color border) {
        setProperties(cursor);
        this.R = R;
        this.B = B;
        this.G = G;
        this.border = border;
    }

    /**
     * Sets the properties of the combo box.
     *
     * @param cursor the cursor to be set for the combo box
     */
    private void setProperties(Cursor cursor) {
        setOpaque(false);
        setCursor(cursor);
    }

    /**
     * Paints the component on the screen.
     *
     * @param g the Graphics object used for painting
     */
    @Override
    protected void paintComponent(Graphics g) {
        g.setColor(new Color(R, B, G));
        g.fillRoundRect(0, 0, getWidth() - 1, getHeight() - 1, getHeight() - 1, getHeight() - 1);
        super.paintComponent(g);
    }

    /**
     * Paints the border of the component.
     *
     * @param g the Graphics object used for painting
     */
    @Override
    protected void paintBorder(Graphics g) {
        g.setColor(border);
        g.drawRoundRect(0, 0, getWidth() - 1, getHeight() - 1, getHeight() - 1, getHeight() - 1);
    }
}
