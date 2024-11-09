from flask_wtf import FlaskForm
from wtforms import StringField, FloatField, SubmitField
from wtforms.validators import DataRequired

class SimulationForm(FlaskForm):
    case = StringField('Simulation Case', default='Custom', validators=[DataRequired()])
    dt = FloatField('Time Step', validators=[DataRequired()])
    H = FloatField('Wave Height', validators=[DataRequired()])
    outputdir = StringField('Output Directory', validators=[DataRequired()])
    submit = SubmitField('Start Simulation')