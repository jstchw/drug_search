import { Container, Row } from 'react-bootstrap';
import Header from '../../components/Header'
import './Settings.css'

const Settings = () => {
    return (
        <div>
            <Header />
            <Container>
                <Row>
                    <h1>Settings</h1>
                </Row>
            </Container>
        </div>
    );
}

export default Settings;