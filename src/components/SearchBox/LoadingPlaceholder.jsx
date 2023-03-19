import { Placeholder } from 'react-bootstrap';

const LoadingPlaceholder = () => {
    return (
        <Placeholder as="div" animation="glow">
            <Placeholder xs={6} />
        </Placeholder>
    );
}

export default LoadingPlaceholder;